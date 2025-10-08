using CSV
using DataFrames
using DelimitedFiles
using Distributions
using LinearAlgebra
using ProgressMeter
using Random
using Statistics;
using Dates;
using JLD2;
using Plots;

data_path = ARGS[1]
folder = ARGS[2]
niter = ARGS[3]
niter = parse(Int64, niter)
@show niter
seed  = ARGS[4]
seed  = parse(Int64, seed)
@show seed

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
@show nInd,nMarker
bv1 = vec(readdlm(data_path * "bv1.txt"));


cluster_size = 3
TransformedX_dict = JLD2.load(data_path * "nrank$cluster_size.eigen/995Eigen/TransformedX_dict.jld2")["TransformedX_dict"]
TransformedY_dict = JLD2.load(data_path * "nrank$cluster_size.eigen/995Eigen/TransformedY_dict.jld2")["TransformedY_dict"]
blkSNPsIndex_dict = JLD2.load(data_path * "nrank$cluster_size.eigen/995Eigen/blkSNPsIndex_dict.jld2")["blkSNPsIndex_dict"]
blkID = Int.(vec(readdlm(data_path * "nrank$cluster_size.eigen/995Eigen/blkIDs.txt", ',')))
nblk = length(blkID)
nsnp = sum(map(length, values(blkSNPsIndex_dict)))
@show nsnp

vary = 1.0
startPi = 0.9

function ST_SBayesCPi_SVD_LDblock(TransformedX_dict, TransformedY_dict; vary, n, nsnp, pi=startPi, niter, estimatePi=true, outFreq=50,seed=seed)
    # R matrix is decomposed to important eigen-values lambda and eigen-vectors
    m = nsnp
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

    Wcorr_dict = deepcopy(TransformedY_dict)
    vare_blk = zeros(nblk)
    varg_blk = zeros(nblk)
    iout = 1

    Random.seed!(seed)
    @showprogress "running MCMC ..." for iter in 1:niter
        nsnp_till_last_blk = 0
        # sampling SNP effects
        logPi = log(pi)
        logPiComp = log(1 - pi)
        invVarEffects = 1 ./ varEffects
        logVarEffects = log.(varEffects)
        logDelta0 = logPi

        for b in 1:nblk
            blk = blkID[b]
            Q = TransformedX_dict[blk]
            wcorr = Wcorr_dict[blk]
            blkSNPsIndex = blkSNPsIndex_dict[blk]
            nsnpb = length(blkSNPsIndex)
            q = length(wcorr)  # number of eigen vectors

            for i in 1:nsnpb
                x = Q[:, i]
                oldAlpha = α[i+nsnp_till_last_blk]
                rhs = (dot(x, wcorr) + oldAlpha) / (vare / n)
                invLhs = 1.0 / (n / vare + invVarEffects)
                uhat = invLhs * rhs
                logDelta1 = 0.5 * (log(invLhs) - logVarEffects + uhat * rhs) + logPiComp
                probDelta1 = 1 / (1 + exp(logDelta0 - logDelta1))

                if (rand() < probDelta1)
                    δ[i+nsnp_till_last_blk] = 1
                    β[i+nsnp_till_last_blk] = uhat + randn() * sqrt(invLhs)
                    α[i+nsnp_till_last_blk] = β[i+nsnp_till_last_blk]
                    wcorr = wcorr + x*(oldAlpha - α[i+nsnp_till_last_blk])
                else
                    if (oldAlpha != 0.0)
                        wcorr = wcorr + x * oldAlpha
                    end
                    δ[i+nsnp_till_last_blk] = 0
                    β[i+nsnp_till_last_blk] = randn() * sqrt(varEffects)
                    α[i+nsnp_till_last_blk] = 0
                end
            end   # end of marker loop 
            nsnp_till_last_blk += nsnpb
            Wcorr_dict[blk] = wcorr
            vare_blk[b] = (n * dot(wcorr, wcorr) + nue * scalee) / rand(Chisq(q + nue))

            # compute genetic variance and heritability
            what = Q * α[blkSNPsIndex]   # w = Qβ + ϵ
            varg_blk[b] = dot(what, what)      # y = Xβ + e 

        end # end of blk loop

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
        vare = mean(vare_blk)
        vare_mcmc[iter] = vare

        # compute genetic variance and heritability 
        varg = sum(varg_blk)     # y = Xβ + e 

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

@time posterior_means = ST_SBayesCPi_SVD_LDblock(TransformedX_dict, TransformedY_dict; vary, n=nInd, nsnp, pi=startPi, niter, estimatePi=true, outFreq=50, seed=seed)

CSV.write(folder * "posterior_means_bcpi.csv", posterior_means)

res_bcpi = posterior_means[7:end, 2]
ebv1_bcpi = geno_qc_df * res_bcpi
accuracy1_bcpi = cor(ebv1_bcpi, bv1)
println("Prediction accuracy:", accuracy1_bcpi)



