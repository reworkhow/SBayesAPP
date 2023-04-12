using CSV
using DataFrames
using DelimitedFiles
using Distributions
using LinearAlgebra
using ProgressMeter
using Random
using Statistics; 
using Dates;

data_path = ARGS[1]
folder = ARGS[2]
niter = ARGS[3]
niter = parse(Int64, niter)
seed = ARGS[4]
seed = parse(Int64, seed)
Random.seed!(seed)
@show seed

#data_path = "/Users/apple/Desktop/JiayiQu/UCD_PhD/MTSBayesC/data/"
#folder = "/Users/apple/Desktop/JiayiQu/UCD_PhD/MTSBayesC/analysis/Sim2T/LDblocks/ST/test/"
#nIter = 1000

# built-in functions 
#single-trait
function samplePi(nEffects::Number, nTotal::Number)
    return rand(Beta(nTotal - nEffects + 1, nEffects + 1))
end
function sample_variance(x, n, df, scale)
    return (dot(x, x) + df * scale) / rand(Chisq(n + df))
end

# Input data 
#read genotypes
geno_qc_df, SNPIDs = readdlm(data_path * "geno_n503_p18588_realSnpID.QC.scale.csv", ',', header=true) # genotypes (after QC)
SNPIDs = SNPIDs[1, 2:end];
geno_qc_df = convert(Matrix{Float64}, geno_qc_df[:, 2:end]);

nInd, nMarker = size(geno_qc_df);
LDmatrix = geno_qc_df'geno_qc_df / (nInd - 1);


### input of function
bhat = vec(readdlm(data_path * "bhat1.txt"));
y1   = vec(readdlm(data_path * "y1.txt"));
bv1  = vec(readdlm(data_path * "bv1.txt"));
### input of function
vary = var(y1)  # variance of y
@show vary
startPi = 0.9
estimatePi = true;

# eigen-decomposition of R matrix
#eigen_values = eigvals(LDmatrix)
#selected = findall(x -> x > 1e-3, eigen_values)
#lambda = eigen_values[selected]
lambda = vec(readdlm(data_path * "eigen_values.txt"));

#eigen_vectors = eigvecs(LDmatrix)
#U = eigen_vectors[:,selected];   # Uq of dimension n_snp x n_eigen
U = readdlm(data_path * "eigen_vectors.txt");

# choose 99.5% of eigenvalues
# Sort the eigenvalues in decreasing order.
decreasing_idx = sortperm(lambda, rev=true)
eigen_values  = lambda[decreasing_idx]
eigen_vectors = U[:, decreasing_idx]
cumulative_sum= cumsum(eigen_values)
# Find the index of the first eigenvalue that is greater than or equal to 99.5% of the total variance.
stop_index        = findfirst(cumulative_sum .>= 0.995 * sum(eigen_values))
lambda_995Eigen   = eigen_values[1:stop_index]
U_995Eigen        = eigen_vectors[:,1:stop_index]


function ST_SBayesCPi_SVD(bhat, U, lambda; vary, n, pi=startPi, niter, estimatePi=true, outFreq=50)
    # R matrix is decomposed to important eigen-values lambda and eigen-vectors
    m = size(U, 1)
    h2 = 0.5
    varg = vary * h2
    vare = vary * (1 - h2)

    varEffects = varg / (m * (1 - pi))
    nub = 4
    nue = 4
    scaleb = (nub - 2) / nub * varEffects
    scalee = (nue - 2) / nue * vare

    α = zeros(m) # a vector of estimated marker effects
    β = zeros(m)
    δ = ones(m)

    α_mcmc = zeros(niter, m)
    pi_mcmc = zeros(niter)
    varEffects_mcmc = zeros(niter)
    vare_mcmc = zeros(niter)
    hsq_mcmc = zeros(niter)

    nOutput = Int(floor(niter / outFreq))
    keptIter = zeros(nOutput, 6)


    w = Diagonal(1 ./ sqrt.(lambda)) * U' * bhat  # w = Qβ + ϵ; var(ϵ) = vare/n I 
    # w = Λq^{-1/2}Uq'b
    Q = Diagonal(sqrt.(lambda)) * U'            # Q = Λq^{1/2}Uq'
    wcorr = deepcopy(w)
    q = size(U, 2)                               # number of eigen vectors 
    iout = 1

    @showprogress "running MCMC ..." for iter in 1:niter
        # sampling SNP effects
        logPi = log(pi)
        logPiComp = log(1 - pi)
        invVarEffects = 1 ./ varEffects
        logVarEffects = log.(varEffects)
        logDelta0 = logPi
        for i in 1:m
            x = Q[:, i]
            oldAlpha = α[i]
            rhs = (dot(x, wcorr) + oldAlpha) / (vare / n)
            invLhs = 1.0 / (n / vare + invVarEffects)
            uhat = invLhs * rhs
            logDelta1 = 0.5 * (log(invLhs) - logVarEffects + uhat * rhs) + logPiComp
            probDelta1 = 1 / (1 + exp(logDelta0 - logDelta1))

            if (rand() < probDelta1)
                δ[i] = 1
                β[i] = uhat + randn() * sqrt(invLhs)
                α[i] = β[i]
                BLAS.axpy!(oldAlpha - α[i], x, wcorr)
            else
                if (oldAlpha != 0.0)
                    BLAS.axpy!(oldAlpha, x, wcorr)
                end
                δ[i] = 0
                β[i] = randn() * sqrt(varEffects)
                α[i] = 0
            end
        end   # end of marker loop 
        α_mcmc[iter, :] = α

        # sampling pi
        nloci = sum(δ)
        if (estimatePi)
            pi = samplePi(nloci, m)
            pi_mcmc[iter] = pi
        end
        # sampling SNP effect variance
        varEffects = sample_variance(α, nloci, nub, scaleb)
        varEffects_mcmc[iter] = varEffects

        # sampling residual variance
        vare = (n * dot(wcorr, wcorr) + nue * scalee) / rand(Chisq(q + nue))
        vare_mcmc[iter] = vare

        # compute genetic variance and heritability
        # ?why don't use vary - vare 
        what = Q * α               # w = Qβ + ϵ
        varg = dot(what, what)      # y = Xβ + e 
        # 2nd way:
        #varg = vary-vare
        #varg_mcmc[iter] = varg

        # heritability 
        hsq = varg / vary
        hsq_mcmc[iter] = hsq

        if iter % outFreq == 0
            keptIter[iout, :] = [pi, nloci, varEffects, hsq, vare, varg]
            println("iter $iter, pi = $pi, nloci = $nloci, varEffects = $varEffects, hsq = $hsq, vare = $vare, varg = $varg.")
            iout += 1
        end
    end # end of iter loop

    df_name = ["pi", "nloci", "varEffects", "hsq", "Vare", "Varg"]
    par = mean(keptIter, dims=1)
    α_mean = mean(α_mcmc, dims=1)

    df = DataFrame(:parameter => vec([df_name; "SNP" .* string.(collect(1:m))]), :mean => vec([par'; α_mean']))

    # save mcmc samples
    writedlm(folder * "mcmcPi.txt", pi_mcmc)
    writedlm(folder * "mcmcVarEffects.txt", varEffects_mcmc)
    writedlm(folder * "mcmcVare.txt", vare_mcmc)
    writedlm(folder * "mcmcHsq.txt", hsq_mcmc)
    writedlm(folder * "mcmc_Alpha.txt", α_mcmc)

    return df
end


@time posterior_means = ST_SBayesCPi_SVD(bhat, U_995Eigen, lambda_995Eigen; vary, n = nInd, pi = startPi, niter)

CSV.write(folder * "posterior_means_bcpi.csv", posterior_means)

res_bcpi = posterior_means[7:end,2]
ebv1_bcpi = geno_qc_df*res_bcpi
accuracy1_bcpi = cor(ebv1_bcpi, bv1)
println("Prediction accuracy:", accuracy1_bcpi)
