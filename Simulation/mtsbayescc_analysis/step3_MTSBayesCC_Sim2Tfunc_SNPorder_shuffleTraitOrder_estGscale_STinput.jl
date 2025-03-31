using CSV
using DataFrames
using DelimitedFiles
using Distributions
using LinearAlgebra
using ProgressMeter
using Random
using Statistics
using Dates
using JLD2
using MPI
using Dates
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--seed"
        help = ""
        arg_type = Int64
        default = 123

        "--pleio_percent"
        help = ""
        arg_type = Float64
        default = 0.1

        "--sample_size"
        help = ""
        arg_type = Int64
        default = 10

        # "--annotation_size"
        # help = ""
        # arg_type = Float64
        # default = 0.1

        "--h21"
        help = ""
        arg_type = Float64
        default = 0.01

        "--h22"
        help = ""
        arg_type = Float64
        default = 0.01

        "--niter"
        help = ""
        arg_type = Int64
        default = 1000

        "--nrank"
        help = ""
        arg_type = Int64
        default = 2

        "--nmarker"
        help = ""
        arg_type = Int64
        default = 100000

        "--outfreq"
        help = ""
        arg_type = Int64
        default = 100

        "--analysis_folder"
        help = ""
        arg_type = String
        default = "./"

        "--ST_folder"
        help = ""
        arg_type = String
        default = "./"

        
    end

    return parse_args(s)
end

# Use the parsed arguments
args = parse_commandline()

##############################################################
## Input data
##############################################################
data_folder = "/common/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/"
## Input data
##############################################################

nInd = args["sample_size"]
# annotation_size = args["annotation_size"]
pleio_percent = args["pleio_percent"]
if pleio_percent == 0 || pleio_percent == 1
    pleio_percent = Int(pleio_percent) #0, instead of 0.0; 1, instead of 1.0
end
seed = args["seed"]
h21 = args["h21"] # h2 for trait 1
h22 = args["h22"] # h2 for trait 2
nIter = args["niter"]
nrank = args["nrank"]
nMarker = args["nmarker"]
outFreq = args["outfreq"]

analysis_folder = args["analysis_folder"]
ST_folder = args["ST_folder"]

#data_folder_name = "h2_trait1.$h21.h2_trait2.$(h22)_pleioPercent$(pleio_percent)_sampleSize$(nInd)_annotationSize$(annotation_size)_seed$seed/"
data_folder_name = "h2_trait1.$h21.h2_trait2.$(h22)_pleioPercent$(pleio_percent)_sampleSize$(nInd)_seed$seed/"
data_path = data_folder * data_folder_name
analysis_path = analysis_folder * data_folder_name
ST_path = ST_folder * data_folder_name
mkpath(analysis_path)
cd(analysis_path)


function sample_variance(ycorr_array, nobs, df, scale, invweights, constraint)
    if invweights != false
        invweights = Diagonal(invweights)
    end
    ntraits = length(ycorr_array)
    SSE = zeros(ntraits, ntraits)
    for traiti = 1:ntraits
        ycorri = ycorr_array[traiti]
        for traitj = traiti:ntraits
            ycorrj = ycorr_array[traitj]
            SSE[traiti, traitj] = (invweights == false) ? dot(ycorri, ycorrj) : ycorri' * invweights * ycorrj
            if constraint == true #diagonal elements only
                break
            end
            SSE[traitj, traiti] = SSE[traiti, traitj]
        end
    end
    if constraint == false
        R = rand(InverseWishart(df + nobs, convert(Array, Symmetric(scale + SSE))))
    else  #diagonal elements only, from scale-inv-⁠χ2
        R = zeros(ntraits, ntraits)
        for traiti = 1:ntraits
            R[traiti, traiti] = (SSE[traiti, traiti] + df * scale[traiti]) / rand(Chisq(nobs + df))
        end
    end
    return R
end

# Flatten the vector of vectors into a sequential vector
function flatten(vec_of_vecs)
    return vcat(map(x -> x[:], vec_of_vecs)...)
end

# Reshape the sequential vector back into a vector of vectors
function unflatten(seq_vec, subvec_length)
    return [seq_vec[i:i+subvec_length-1] for i in 1:subvec_length:length(seq_vec)]
end

# Flatten the vector of matrices into a sequential vector
function flatten_matrices(vec_of_mats)
    return vcat(map(x -> x[:], vec_of_mats)...)
end

# Reshape the sequential vector back into a vector of matrices
function unflatten_matrices(seq_vec, nrows, ncols)
    num_matrices = length(seq_vec) ÷ (nrows * ncols)
    return [reshape(seq_vec[(i-1)*nrows*ncols+1:i*nrows*ncols], nrows, ncols) for i in 1:num_matrices]
end

function compute_correlation(covMatrix)
    # Compute the standard deviations (sqrt of diagonal elements)
    std_devs = sqrt.(diag(covMatrix))
    # Compute the correlation coefficient
    correlation = covMatrix[1, 2] / (std_devs[1] * std_devs[2])
    return correlation
end

############################################################################################################
# Input data 
nBlocks = 49
annotationName = ["category1", "category2"]; # ordered by continuous category then categorical category

# SNPc1 = readdlm(data_path * "SNPc1.txt",)[:, 1]
SNPc1 = readdlm(data_folder * "SNPc1.txt",)[:, 1]
nSNPc1 = length(SNPc1)
# SNPc2 = readdlm(data_path * "SNPc2.txt",)[:, 1]
SNPc2 = readdlm(data_folder * "SNPc2.txt",)[:, 1]
nSNPc2 = length(SNPc2)
nLoci_annot = [nSNPc1, nSNPc2]

nCon = 0 # number of continuous annotation
nCat = length(annotationName) # number of categorical annotation
annotationType = repeat(["category"], nCat) # ordered by "continue" then "category"
estimate_vare = true
estimate_vara = true
estimate_pi = true
sample_delta = false
estimate_Gscale = true
############################################################################################################

R = [1 0.1; 0.1 1] / nInd                             # R_star, residual variance

# Pct00 = 0.99 #always 1% QTLs
# Pct11 = (1 - Pct00) * pleio_percent
# Pct10 = Pct01 = (1 - Pct00 - Pct11) * 0.5
#startPi = Dict([1.0; 1.0] => Pct11, [1.0; 0.0] => Pct10, [0.0; 1.0] => Pct01, [0.0; 0.0] => Pct00)
#if pleio_percent == 1 || pleio_percent == 0
#     println("!!! pleio_percent=0 or 100%!!!, reset Pi to avoid 0 values in Pi")
#     startPi = Dict([1.0; 1.0] => 0.01, [1.0; 0.0] => 0.01, [0.0; 1.0] => 0.01, [0.0; 0.0] => 0.97)
# end

#marker effect variance-covariance matrix
# Gprior_vec = [zeros(2, 2) for c in 1:nCat]
# saved_sdbv = readdlm(data_path * "sdbv.txt")[:, 1]
# saved_sdy = readdlm(data_path * "sdy.txt")[:, 1]

# for c in 1:nCat
#     bv_var_file_path = data_path * "bvc$(c)_var_cov.txt"
#     saved_bv_var = readdlm(bv_var_file_path)
#     t1_scaler = 1 / saved_sdbv[1] * sqrt(h21) / saved_sdy[1]
#     t2_scaler = 1 / saved_sdbv[2] * sqrt(h22) / saved_sdy[2]
#     Gprior_vec[c][1,1] = saved_bv_var[1,1] * t1_scaler^2
#     Gprior_vec[c][2,2] = saved_bv_var[2,2] * t2_scaler^2
#     Gprior_vec[c][1,2] = Gprior_vec[c][2,1] = saved_bv_var[1,2] * t1_scaler * t2_scaler
#     Gprior_vec[c] = Gprior_vec[c] / (nLoci_annot[c] * (1 - Pct00))
# end

# for c in 1:nCat
#     marker_var_file_path = data_path * "QTLcovmatc$(c).txt"
#     if isfile(marker_var_file_path)
#         #Gprior_vec[c] = reshape(readdlm(marker_var_file_path), 2, 2)
#         Gprior_vec[c] = readdlm(marker_var_file_path)
#     else
#         error("The file $marker_var_file_path does not exist!!!")
#     end
# end

# estGscale
Pi11 = 0.001
# Read and round values for Trait1 and Trait2
Pi10 = 1.0 - round(readdlm(ST_path * "Trait1/mean_pi.txt")[1,1], digits=3)
Pi01 = 1.0 - round(readdlm(ST_path * "Trait2/mean_pi.txt")[1,1], digits=3)

# Calculate Pi00
Pi00 = 1.0 - Pi11 - Pi10 - Pi01
startPi = Dict([1.0; 1.0] => Pi11, [1.0; 0.0] => Pi10, [0.0; 1.0] => Pi01, [0.0; 0.0] => Pi00)

ST_h21 = round(readdlm(ST_path * "Trait1/mean_varg_total.txt")[1,1], digits=2)
ST_h22 = round(readdlm(ST_path * "Trait2/mean_varg_total.txt")[1,1], digits=2)

Gprior_vec = [zeros(2, 2) for c in 1:nCat]
# for c in 1:nCat
#     Gprior_vec[c] = [ST_h21 0.0; 0.0 ST_h22] / nCat
#     Gprior_vec[c] = Gprior_vec[c] / (nLoci_annot[c] * (1 - Pi00))
# end

for c in 1:nCat
    Gprior_vec[c] = [ST_h21 0.0; 0.0 ST_h22] * (nLoci_annot[c] / nMarker)
    Gprior_vec[c] = Gprior_vec[c] / (nLoci_annot[c] * (1 - Pi00))
end

function runMPI(; 
    startPi=startPi,
    nIter=nIter, outFreq=outFreq, seed=seed, burnin=Int(0.4 * nIter), 
    estimate_vare=estimate_vare, estimate_vara=estimate_vara,
    nInd=nInd, nMarker=nMarker, nBlocks=nBlocks,
    annotationType=annotationType, annotationName=annotationName,
    nCon=nCon, nCat=nCat,
    analysis_path=analysis_path, data_path=data_path,
    Gprior_vec=Gprior_vec, R=R)

    comm = MPI.COMM_WORLD
    my_rank = MPI.Comm_rank(comm) #current rank, e.g., 0/1/2/3, root=0
    cluster_size = MPI.Comm_size(comm) #number of all processes/rank, e.g., 4
    @show my_rank, cluster_size
    MPI.Barrier(comm)

    # set seed in different rank
    Random.seed!(seed + my_rank)
    MPI.Barrier(comm)

    ############################################################################
    #read data in current rank
    ############################################################################

    my_TransformedX_dict = JLD2.load(data_path * "nrank$cluster_size.eigen/bhatXsj/995Eigen/rank$my_rank.TransformedX_dict.jld2")["my_TransformedX_dict"]
    my_TransformedY_dict = JLD2.load(data_path * "nrank$cluster_size.eigen/bhatXsj/995Eigen/rank$my_rank.TransformedY_dict.jld2")["my_TransformedY_dict"] # blk => [y1,y2]
    my_blkSNPsIndex_dict = JLD2.load(data_path * "nrank$cluster_size.eigen/bhatXsj/995Eigen/rank$my_rank.blkSNPsIndex_dict.jld2")["my_blkSNPsIndex_dict"]
    my_blkID = Int.(vec(readdlm(data_path * "nrank$cluster_size.eigen/bhatXsj/995Eigen/rank$my_rank.blkIDs.txt", ',')))  #block IDs for this rank
    my_nblk = length(my_blkID)
    my_nsnp = sum(map(length, values(my_blkSNPsIndex_dict)))
    # need to change accordingly
    my_anno_matrix_dict = JLD2.load(data_path * "nrank$cluster_size.eigen/bhatXsj/995Eigen/rank$my_rank.anno_matrix_dict.jld2")["my_anno_matrix_dict"]

    
    # hyper-parameters for A (marker effect variance) and R
    nCategory = nCat + nCon
    A_vec = Gprior_vec # marker effect variance matrix
    Ainv_vec = [inv(A_vec[c]) for c in 1:nCategory]
    Pi = [copy(startPi) for c in 1:nCategory]

    nTraits = 2
    if estimate_vare == true
        df_R = 4 + nTraits
        scale_R = R * (df_R - nTraits - 1)
    end

    if my_rank == 0
        if estimate_vara == true
            df_G = 4 + nTraits
            scale_G_vec = [Gprior_vec[c] * (df_G - nTraits - 1) for c in 1:nCategory]
        end
    end

    if my_rank == 0
        writedlm(analysis_path * "annotationName.txt", annotationName)
    end
    # reorder my_blkSNPsIndex_dict as from 1 to my_nsnp, before was from 1 to nsnp_per_block
    sort!(my_blkID) # sort my_blkID
    for i in 2:my_nblk
        previousb = my_blkID[i-1]
        ncumSNP = length(my_blkSNPsIndex_dict[previousb])
        for key in my_blkID[i:my_nblk]
            my_blkSNPsIndex_dict[key] = my_blkSNPsIndex_dict[key] .+ ncumSNP
        end
    end

    xpx_dict = Dict{Int,Vector{Vector{Float64}}}() # blkID -> vector of xpx for each category
    xArray_dict = Dict() # blkID -> TransformedX for each category

    for blk in keys(my_TransformedX_dict)
        Xb = my_TransformedX_dict[blk]
        nMarkerb = size(Xb, 2)
        # save xpx for continue and category annotation
        # the last element is for category annotation
        xpx_dict[blk] = Vector{Vector{Float64}}(undef, nCon + 1)
        xArray_dict[blk] = Vector{Matrix{Float64}}(undef, nCon + 1)
        for c in 1:nCon # xpx and xArray will change for continuous group 
            xpx_dict[blk][c] = [my_anno_matrix_dict[blk][:, c][i]^2 * dot(Xb[:, i], Xb[:, i]) for i in 1:nMarkerb]
            xArray_dict[blk][c] = Xb * diagm(my_anno_matrix_dict[blk][:, c])
        end
        # xpx/xArray will not change for categorical group if the SNPs are not in the group (marker effects = 0)
        xpx_dict[blk][end] = [dot(Xb[:, i], Xb[:, i]) for i in 1:nMarkerb]
        xArray_dict[blk][end] = Xb
    end

    # change my_anno_matrix_dict to a true/false matrix dict to indicate whether SNPs to be sampled for each category or not 
    my_anno_matrix_dict = Dict(key => (my_anno_matrix_dict[key] .!= 0.0) for key in keys(my_anno_matrix_dict))
    if nCon > 0
        # for continuous annotation, all values = true
        concol = findall(annotationType .== "continue") # column index for continuous annotation
        for blk in keys(my_anno_matrix_dict)
            my_anno_matrix_dict[blk][:, concol] .= true
        end
    end

    #output
    β = zeros(nTraits)
    newα = zeros(nTraits)  #α=Dβ
    oldα = zeros(nTraits)
    w = zeros(nTraits)
    δ = zeros(nTraits)

    betaArray = [zeros(my_nsnp * nCategory) for t in 1:nTraits] #-> ordered by annotation groups 
    alphaArray = [zeros(my_nsnp * nCategory) for t in 1:nTraits]
    deltaArray = [zeros(my_nsnp * nCategory) for t in 1:nTraits]

    

    #posterior mean 
    meanAlpha = [zeros(my_nsnp * nCategory) for t in 1:nTraits]
    if sample_delta
        mcmc_Delta = [zeros(my_nsnp * nCategory, nIter - burnin) for t in 1:nTraits]
    else
        nOutput = Int(floor((nIter-burnin) / outFreq))
        mcmc_Delta = [zeros(my_nsnp * nCategory, nOutput) for t in 1:nTraits]
    end

    if my_rank == 0
        if estimate_pi == true
            mean_pi = deepcopy(Pi)
            mean_pi2 = deepcopy(Pi)
            for i in 1:nCategory
                for key in keys(mean_pi[i])
                    mean_pi[i][key] = 0.0
                    mean_pi2[i][key] = 0.0
                end
            end
        end

        if estimate_vara == true

            # beta effect variance
            meanB2 = [zeros(nTraits, nTraits) for c in 1:nCategory]
            meanBcor2 = zeros(nCategory)
        end

        #dot(\alpha, \alpha)
        meanSSQ = [zeros(nTraits, nTraits) for c in 1:nCategory]
        meanSSQ2 = [zeros(nTraits, nTraits) for c in 1:nCategory]

        # marker effect variance
        meanA = [zeros(nTraits, nTraits) for c in 1:nCategory]
        meanA2 = [zeros(nTraits, nTraits) for c in 1:nCategory]
        meanAcor = zeros(nCategory)
        meanAcor2 = zeros(nCategory)

        meanB = [zeros(nTraits, nTraits) for c in 1:nCategory]
        meanBcor = zeros(nCategory)

        # genetic effect variance
        meanG = [zeros(nTraits, nTraits) for c in 1:nCategory]
        meanG2 = [zeros(nTraits, nTraits) for c in 1:nCategory]
        meanGcor = zeros(nCategory)
        meanGcor2 = zeros(nCategory)

        meanGenr_ssq = zeros(nCategory)
        meanGenr2_ssq = zeros(nCategory)
        meanGenr_whq = zeros(nCategory)
        meanGenr2_whq = zeros(nCategory)

        meanGtotal = zeros(nTraits, nTraits)
        meanGtotal2 = zeros(nTraits, nTraits)

        mcmcAtruecor_c = zeros(nIter - burnin, nCategory)
        mcmcBcor_c = zeros(nIter - burnin, nCategory)
        mcmcGcov_c = zeros(nIter - burnin, nCategory)
        mcmcGcor_c = zeros(nIter - burnin, nCategory)
        mcmcGenr_ssq_c = zeros(nIter - burnin, nCategory)
        mcmcGenr_whq_c = zeros(nIter - burnin, nCategory)

        mcmcGcov_total = zeros(nIter - burnin)
        mcmcGcor_total = zeros(nIter - burnin)

        if estimate_vare == true
            meanR = zeros(nTraits, nTraits)
            meanR2 = zeros(nTraits, nTraits)
        end
    end

    if my_rank == 0
        # outfile to save mcmc results 
        outfile = Dict{String,IOStream}()
        outvar = ["pi"]
        #if estimate_vara == true
        push!(outvar, "marker_effects_variance")
        push!(outvar, "beta_effects_variance")
        push!(outvar, "genetic_effects_variance")
        #end
        file_name = "MCMC_samples"
        for i in outvar
            file_i = analysis_path * file_name * "_$i.txt"
            if isfile(file_i)
                printstyled("The file " * file_i * " already exists!!! It is overwritten by the new output.\n", bold=false, color=:red)
            else
                printstyled("The file " * file_i * " is created to save MCMC samples for " * i * ".\n", bold=false, color=:green)
            end
            outfile[i] = open(file_i, "w")
        end
    end

    if my_rank == 0
        println("---------------- Summary Start --------------")
        println("nIter=$nIter, outFreq=$outFreq, seed=$seed, burnin = $burnin")
        #println("Annotation size: ", annotation_size)
        println("Pleiotropy percent: ", pleio_percent)
        println("Heritability for trait 1 (h21): ", ST_h21)
        println("Heritability for trait 2 (h22): ", ST_h22)
        println("startPi=$startPi")
        println("Number of ranks: ", nrank)
        println("estimate_vare=$estimate_vare,estimate_vara=$estimate_vara")
        println("estimate_pi=$estimate_pi")
        println("estimate_Gscale=$estimate_Gscale")
        println("analysis_path=$analysis_path")
        println("data_path=$data_path")
        println("marker effect variance for C1 is: ", Gprior_vec[1])
        println("marker effect variance for C2 is: ", Gprior_vec[2])
        println("residual variance is: ", round.(R * nInd, digits=3), "  (R*Ind)")
        println("nMarker=$nMarker, nInd=$nInd, nBlocks=$nBlocks")
        println("nCat = $nCat, nCon = $nCon")
        time_start = now()
        println("Start time: ", time_start)
        println("---------------- Summary End ----------------")
    end
    MPI.Barrier(comm)
    println("In rank$my_rank, there are $my_nblk LD blocks, and $my_nsnp SNPs in total.")
    MPI.Barrier(comm)

    iout = 1
    R_blk = [R for _ in 1:my_nblk]


    @showprogress "running MCMC ..." for iter = 1:nIter
        if iter > burnin
            iIter = 1.0 / (iter - burnin)
        end
        if estimate_Gscale
            # iIter_scaleG is used to sample G_scale (meanGtotal & mean_pi)
            if iter <= burnin
                iIter_scaleG = 1.0 / iter
            end
        end

        # varg Computed as dot(X\alpha, X\alpha)         
        varg_blk_cat = [zeros(nTraits, nCategory) for _ in 1:my_nblk] # varg for different category & trait (saved for hsq computation)
        varg_cov_blk_cat = [zeros(nCategory) for _ in 1:my_nblk] # genetic covariance for different category
        # compute total varg for each trait without split into different categories
        totalvarg_blk = [zeros(nTraits, nTraits) for _ in 1:my_nblk]

        # ssq_blk_cat Computed as dot(\alpha, \alpha)
        ssq_blk_cat = [zeros(nTraits, nCategory) for _ in 1:my_nblk] # sum of square for different category & trait (saved for enrichment)
        ssq_cov_blk_cat = [zeros(nCategory) for _ in 1:my_nblk] # sum of square for different category (saved for enrichment)

        nLoci_array_vec = [fill(0, length(startPi)) for c in 1:nCategory] # number of loci (pi) in each category for SNPs in this rank
        # SSE_vec Computed as dot(\beta, \beta) for whole rank
        SSE_vec = [zeros(nTraits, nTraits) for c in 1:nCategory]
        # G_vec Computed as dot(X\alpha, X\alpha) for whole rank
        G_vec = [zeros(nTraits, nTraits) for c in 1:nCategory] # genetic variance for different category for SNPs in this rank

        
        #for loop for LD blocks
        for b in 1:my_nblk
            blk = my_blkID[b] #block ID
            xpx_vec = xpx_dict[blk]
            xArray_vec = xArray_dict[blk]
            wArray = my_TransformedY_dict[blk]
            annotationMatb = my_anno_matrix_dict[blk]  #annotation boolean data for current blk, nMarkerb-by-C
            SNPIndexb = my_blkSNPsIndex_dict[blk]
            nEigenb, nMarkerb = size(my_TransformedX_dict[blk]) # q & nsnpb

            # # initialize xArrayc/xpxc 
            xArrayc = xArray_vec[1] # nEigenb x nMarkerb matrix
            xpxc = xpx_vec[1]

            R = R_blk[b]
            Rinv = inv(R)

            
            for marker = 1:nMarkerb
                true_marker_num = SNPIndexb[marker] #marker position across my_nsnp
                annoindexm = annotationMatb[marker, :] # categories for current marker

                for cat = 1:nCategory
                    if nCon != 0
                        if cat <= nCon + 1 # continuous group + categorical group (after nCon+1, xArrayc/xpxc will not change)
                            xArrayc = xArray_vec[cat]
                            xpxc = xpx_vec[cat]
                        end
                    end
                    Ginv = Ainv_vec[cat]
                    BigPi = Pi[cat]

                    if annoindexm[cat] # skip sampling if the marker is not in the category
                        markerIndex = (cat - 1) * my_nsnp + true_marker_num # exact position in betaArray; betaArray[trait1]: nMarker*nCategory-by-1
                        x = xArrayc[:, marker]
                        for trait = 1:nTraits
                            β[trait] = betaArray[trait][markerIndex]
                            oldα[trait] = newα[trait] = alphaArray[trait][markerIndex]
                            δ[trait] = deltaArray[trait][markerIndex]
                            w[trait] = dot(x, wArray[trait]) + xpxc[marker] * oldα[trait] #w=xj'(ycorr+xj*αj) scaler
                        end


                        TraitOrder = shuffle(1:nTraits)
                        for k = TraitOrder
                            Ginv11 = Ginv[k, k]
                            nok = deleteat!(collect(1:nTraits), k)
                            Ginv12 = Ginv[k, nok] #this is not row vector!!, so this is Ginv21
                            C11 = Ginv11 + Rinv[k, k] * xpxc[marker]
                            C12 = Ginv12 + xpxc[marker] * Diagonal(δ[nok]) * Rinv[k, nok] #C21
                            #when δj=0
                            invLhs0 = 1 / Ginv11
                            rhs0 = -dot(Ginv12, β[nok])
                            gHat0 = rhs0 * invLhs0
                            #when δj=1
                            invLhs1 = 1 / C11
                            rhs1 = w' * Rinv[:, k] - C12'β[nok]  #here the w' is in paper: w'm=xj'(ycorr+xj*βj)
                            gHat1 = rhs1 * invLhs1

                            d0 = copy(δ)
                            d1 = copy(δ)
                            d0[k] = 0.0
                            d1[k] = 1.0

                            #sample δj
                            logDelta0 = -0.5 * (log(Ginv11) - gHat0^2 * Ginv11) + log(BigPi[d0]) #logPi
                            logDelta1 = -0.5 * (log(C11) - gHat1^2 * C11) + log(BigPi[d1]) #logPiComp
                            probDelta1 = 1.0 / (1.0 + exp(logDelta0 - logDelta1))

                            #sample marker effects
                            if (rand() < probDelta1) # force all goes to [1,1]
                                δ[k] = 1
                                β[k] = newα[k] = gHat1 + randn() * sqrt(invLhs1)
                                wArray[k] = wArray[k] + x * (oldα[k] - newα[k])
                            else
                                β[k] = gHat0 + randn() * sqrt(invLhs0)
                                δ[k] = 0
                                newα[k] = 0
                                if oldα[k] != 0
                                    wArray[k] = wArray[k] + x * oldα[k] #newα[k]=0
                                end
                            end
                        end

                        # add to nLoci_array_vec based on δ
                        pi_index = 1
                        for key in keys(BigPi)
                            if δ == key
                                nLoci_array_vec[cat][pi_index] += 1
                            end
                            pi_index += 1
                        end

                        for trait = 1:nTraits
                            betaArray[trait][markerIndex] = β[trait]
                            deltaArray[trait][markerIndex] = δ[trait]
                            alphaArray[trait][markerIndex] = newα[trait]
                            if iter > burnin
                                meanAlpha[trait][markerIndex] += (newα[trait] - meanAlpha[trait][markerIndex]) * iIter
                            end
                        end
                    end # end if loop  
                end # end annotation loop
            end  # end marker loop

            my_TransformedY_dict[blk] = wArray
            
            # compute genetic variance and heritability 
            XAb = xArray_vec[1]
            what_array_total = [zeros(size(XAb, 1)) for i in 1:nTraits] # used to compute total genetic variance
            
            for cat in 1:nCategory
                if nCon != 0
                    if cat <= nCon + 1 # continuous group + categorical group (after nCon+1, xArrayc/xpxc will not change)
                        XAb = xArray_vec[cat]
                    end
                end
                what_array_c = [zeros(size(XAb, 1)) for i in 1:nTraits] # used to compute genetic variance for different category & trait
                alphaArray_c = [alphaArray[i][(cat-1)*my_nsnp.+SNPIndexb] for i in 1:nTraits]
                for traiti = 1:nTraits
                    what_array_c[traiti][:] = XAb * alphaArray_c[traiti]
                    varg_blk_cat[b][traiti, cat] = dot(what_array_c[traiti], what_array_c[traiti]) # genetic variance for different category & trait in bth block
                    ssq_blk_cat[b][traiti, cat] = dot(alphaArray_c[traiti], alphaArray_c[traiti])

                    what_array_total[traiti] += what_array_c[traiti]  # add up whati across category to compute total genetic variance
                end
                varg_cov_blk_cat[b][cat] = dot(what_array_c[1], what_array_c[2]) # genetic covariance for different category in bth block
                ssq_cov_blk_cat[b][cat] = dot(alphaArray_c[1], alphaArray_c[2])
            end
        

            # total genetic variance
            for traiti in 1:nTraits
                for traitj in traiti:nTraits
                    totalvarg_blk[b][traiti, traitj] = dot(what_array_total[traiti], what_array_total[traitj])
                    totalvarg_blk[b][traitj, traiti] = totalvarg_blk[b][traiti, traitj]
                end
            end
        

            if estimate_vare == true
                # sample bivariate residual variance R 
                R_blk[b] = sample_variance(wArray, nEigenb, df_R, scale_R, false, false)
            end
        end # end block loop 

        # summing residual variance
        if estimate_vare == true
            R_blk_sum_rank = sum(R_blk)
        end

        # get \alpha'\alpha to compute enrichment 
        if iter > burnin
            SSQ_vec = [zeros(nTraits, nTraits) for c in 1:nCategory]
            ssq_cat = sum(ssq_blk_cat)
            ssq_cov_cat = sum(ssq_cov_blk_cat)

            for cat = 1:nCategory
                SSQ_vec[cat][1, 1] = ssq_cat[1, cat]
                SSQ_vec[cat][2, 2] = ssq_cat[2, cat]
                SSQ_vec[cat][1, 2] = SSQ_vec[cat][2, 1] = ssq_cov_cat[cat]
            end

            SSQ_vec_flatten = flatten_matrices(SSQ_vec)
            SSQ_vec_flatten_sum = MPI.Reduce(SSQ_vec_flatten, +, 0, comm)
            MPI.Barrier(comm)

            if my_rank == 0
                SSQ_vec_sum = unflatten_matrices(SSQ_vec_flatten_sum, nTraits, nTraits)
                for cat = 1:nCategory
                    meanSSQ[cat] += (SSQ_vec_sum[cat] - meanSSQ[cat]) * iIter
                    meanSSQ2[cat] += (SSQ_vec_sum[cat] .^ 2 - meanSSQ2[cat]) * iIter
                end
            end
        end
        # get true A matrix by alphaArray
        # use only pleiotropic markers to compute QTL effect variance matrix
        if iter > burnin
            Atrue_vec = [zeros(nTraits, nTraits) for c in 1:nCategory]
            nQTL = zeros(nCategory)
            for cat = 1:nCategory
                alpha_array = [alphaArray[i][((cat-1)*my_nsnp+1):((cat-1)*my_nsnp+my_nsnp)] for i in 1:nTraits]
                pleio_marker = (alpha_array[1] .!= 0) .& (alpha_array[2] .!= 0)
                nQTL[cat] = sum(pleio_marker)
                for traiti in 1:nTraits
                    for traitj in traiti:nTraits
                        Atrue_vec[cat][traiti, traitj] = dot(alpha_array[traiti][pleio_marker], alpha_array[traitj][pleio_marker])
                        Atrue_vec[cat][traitj, traiti] = Atrue_vec[cat][traiti, traitj]
                    end
                end
            end
            Atrue_vec_flatten = flatten_matrices(Atrue_vec)
            # gather the results from all ranks
            Atrue_vec_flatten_sum = MPI.Reduce(Atrue_vec_flatten, +, 0, comm)
            nQTL_sum = MPI.Reduce(nQTL, +, 0, comm)
            if my_rank == 0
                Atrue_vec_sum = unflatten_matrices(Atrue_vec_flatten_sum, nTraits, nTraits)
                for cat = 1:nCategory
                    #annotation specific A
                    if nQTL_sum[cat] == 0
                        Atrue_cat = zeros(nTraits, nTraits)
                    else
                        Atrue_cat = Atrue_vec_sum[cat] / nQTL_sum[cat]
                    end
                    meanA[cat] += (Atrue_cat - meanA[cat]) * iIter
                    meanA2[cat] += (Atrue_cat .^ 2 - meanA2[cat]) * iIter
                    # correlation values
                    mcmcAtruecor_c[iter-burnin, cat] = compute_correlation(Atrue_cat)
                end
            end
        end
        MPI.Barrier(comm)
        ########################
        ### Step2. sample Pi ###
        ########################
        # gather nLoci_array_vec from all ranks 
        # flatten nLoci_array_vec
        nLoci_array_vec_flatten = flatten(nLoci_array_vec)
        nLoci_array_vec_flatten_sum = MPI.Reduce(nLoci_array_vec_flatten, +, 0, comm) # sum the results from all ranks
        if my_rank == 0
            nLoci_array_vec_sum = unflatten(nLoci_array_vec_flatten_sum, length(startPi))
        end
        if estimate_pi == true
            tempPi_vec = [zeros(length(startPi)) for cat = 1:nCategory]
            if my_rank == 0
                for cat = 1:nCategory
                    tempPi_vec[cat] = rand(Dirichlet(nLoci_array_vec_sum[cat] .+ 1))
                    if estimate_Gscale || iter > burnin
                        iIter_Pi = estimate_Gscale && iter <= burnin ? iIter_scaleG : iIter
                        tempPi2 = tempPi_vec[cat] .^ 2
                        iCategori = 1
                        for i in keys(Pi[cat])
                            mean_pi[cat][i] += (tempPi_vec[cat][iCategori] - mean_pi[cat][i]) * iIter_Pi
                            mean_pi2[cat][i] += (tempPi2[iCategori] - mean_pi2[cat][i]) * iIter_Pi
                            iCategori += 1
                        end
                    end
                end
            end
        
            ############################################################
            # broadcast tempPi_vec from 0 to other ranks 
            tempPi_vec_flatten = flatten(tempPi_vec)
            tempPi_vec_flatten = MPI.bcast(tempPi_vec_flatten, 0, comm)
            tempPi_vec = unflatten(tempPi_vec_flatten, length(startPi))
            for cat = 1:nCategory
                iCategori = 1
                for i in keys(Pi[cat])
                    Pi[cat][i] = tempPi_vec[cat][iCategori] #annotation specific pi
                    iCategori = iCategori + 1
                end
            end
        end
        MPI.Barrier(comm)

        # summing genetic variance
        if estimate_Gscale || iter > burnin
            varg_cat = sum(varg_blk_cat) # sum varg_blk_cat to get varg for different category across blocks in this rank
            varg_cov_cat = sum(varg_cov_blk_cat) # sum varg_cov_blk_cat to get genetic covariance for different category across blocks in this rank

            for cat = 1:nCategory
                G_vec[cat][1, 1] = varg_cat[1, cat]
                G_vec[cat][2, 2] = varg_cat[2, cat]
                G_vec[cat][1, 2] = G_vec[cat][2, 1] = varg_cov_cat[cat]
            end

            G_vec_flatten = flatten_matrices(G_vec)
            # sum the G_vec_flatten from all ranks to root rank
            G_vec_flatten_sum = MPI.Reduce(G_vec_flatten, +, 0, comm)
            MPI.Barrier(comm)

            if my_rank == 0
                G_vec_sum = unflatten_matrices(G_vec_flatten_sum, nTraits, nTraits)

                # Determine iIter_G based on conditions
                iIter_G = if estimate_Gscale && iter <= burnin
                    iIter_scaleG
                else
                    iIter
                end

                # Calculate meanG and meanG2 regardless of estimate_Gscale
                for cat = 1:nCategory
                    meanG[cat] += (G_vec_sum[cat] - meanG[cat]) * iIter_G
                    meanG2[cat] += (G_vec_sum[cat] .^ 2 - meanG2[cat]) * iIter_G
                end              
            end

            if iter > burnin
                totalvarg = sum(totalvarg_blk) # sum totalvarg_blk to get total genetic variance across all blocks in this rank
                G_total = MPI.Reduce(totalvarg, +, 0, comm)
                MPI.Barrier(comm)
                
                if my_rank == 0
                    mcmcGcor_total[iter-burnin] = compute_correlation(G_total)
                    mcmcGcov_total[iter-burnin] = G_total[1, 2]
                    SSQ_total = sum(SSQ_vec_sum)
                    
                    for cat = 1:nCategory
                        # correlation values
                        mcmcGcor_c[iter-burnin, cat] = compute_correlation(G_vec_sum[cat])
                        mcmcGcov_c[iter-burnin, cat] = G_vec_sum[cat][1, 2]
                        # two ways to compute enrichment
                        mcmcGenr_ssq_c[iter-burnin, cat] = SSQ_vec_sum[cat][1, 2] / SSQ_total[1, 2]
                        mcmcGenr_whq_c[iter-burnin, cat] = G_vec_sum[cat][1, 2] / G_total[1, 2]
                    end
                    meanGtotal += (G_total - meanGtotal) * iIter
                    meanGtotal2 += (G_total .^ 2 - meanGtotal2) * iIter
                end
            end
        end
        MPI.Barrier(comm)
        ############################################################
        ########################
        ### sample scale_G ###
        ########################
        ########################
        ### sample scale_G ###
        ########################
        ########################
        ### sample scale_G ###
        ########################
        if estimate_Gscale
            if iter <= burnin
                if my_rank == 0
                    # meanG
                    # mean_pi
                    Gprior_vec = copy(meanG)
                    for cat = 1:nCategory
                        Gprior_vec[cat][1,1] = meanG[cat][1,1]/(nLoci_annot[cat] * (mean_pi[cat][[1.0, 1.0]] + mean_pi[cat][[1.0, 0.0]]))
                        Gprior_vec[cat][2,2] = meanG[cat][2,2]/(nLoci_annot[cat] * (mean_pi[cat][[1.0, 1.0]] + mean_pi[cat][[0.0, 1.0]]))
                        Gprior_vec[cat][1,2] = Gprior_vec[cat][2,1] = meanG[cat][1,2]/(nLoci_annot[cat] * (mean_pi[cat][[1.0, 1.0]]))
                        scale_G_vec[cat] = Gprior_vec[cat] * (df_G - nTraits - 1)
                    end

                    if iter == burnin
                        # save the scale_G_vec
                        for cat = 1:nCategory
                            writedlm(analysis_path * "scale_G" * string(cat) * ".txt", scale_G_vec[cat])
                        end
                        # re-initialize meanG and mean_pi
                        meanG = [zeros(nTraits, nTraits) for c in 1:nCategory]
                        meanG2 = [zeros(nTraits, nTraits) for c in 1:nCategory]
                        if estimate_pi
                            mean_pi = deepcopy(Pi)
                            mean_pi2 = deepcopy(Pi)
                            for i in 1:nCategory
                                for key in keys(mean_pi[i])
                                    mean_pi[i][key] = 0.0
                                    mean_pi2[i][key] = 0.0
                                end
                            end
                        end
                    end
                end   
            end
        end

        ### Step3. sample variance ###
        # get SSE_vec to sample A in root
        # assume estimate_vara = true
        # add to SSE_vec based on β (marker effect variance)
        for cat in 1:nCategory
            for traiti in 1:nTraits
                beta_i = betaArray[traiti][((cat-1)*my_nsnp+1):((cat-1)*my_nsnp+my_nsnp)]
                for traitj in traiti:nTraits
                    beta_j = betaArray[traitj][((cat-1)*my_nsnp+1):((cat-1)*my_nsnp+my_nsnp)]
                    SSE_vec[cat][traiti, traitj] = dot(beta_i, beta_j)
                    SSE_vec[cat][traitj, traiti] = SSE_vec[cat][traiti, traitj]
                end
            end
        end
        SSE_vec_flatten = flatten_matrices(SSE_vec)

        MPI.Barrier(comm) 
        # gather the results from all ranks
        SSE_vec_flatten_sum = MPI.Reduce(SSE_vec_flatten, +, 0, comm)
        if my_rank == 0
            SSE_vec_sum = unflatten_matrices(SSE_vec_flatten_sum, nTraits, nTraits)
            for cat = 1:nCategory
                if estimate_vara == true
                    A_vec[cat] = rand(InverseWishart(df_G + sum(nLoci_array_vec_sum[cat]), convert(Array, Symmetric(scale_G_vec[cat] + SSE_vec_sum[cat]))))
                end
                if iter > burnin
                    # save mean for beta effect variance 
                    meanB[cat] += (A_vec[cat] - meanB[cat]) * iIter
                    # correlation values
                    mcmcBcor_c[iter-burnin, cat] = compute_correlation(A_vec[cat])
                    if estimate_vara == true
                        meanB2[cat] += (A_vec[cat] .^ 2 - meanB2[cat]) * iIter
                    end
                end
            end
        end
        if estimate_vara == true  
            # broadcast A_vec from 0 to other ranks
            A_vec_flatten = flatten_matrices(A_vec)
            A_vec_flatten = MPI.bcast(A_vec_flatten, 0, comm)
            A_vec = unflatten_matrices(A_vec_flatten, nTraits, nTraits)
            Ainv_vec = [inv(A_vec[cat]) for cat in 1:nCategory]
            MPI.Barrier(comm)
        end

        

        if estimate_vare == true
            # sampling residual variance
            R_blk_sum = MPI.Reduce(R_blk_sum_rank, +, 0, comm)
            MPI.Barrier(comm)
            if my_rank == 0
                R_blkmean = R_blk_sum / nBlocks
                if iter > burnin
                    R2 = (R_blkmean * nInd) .^ 2 #blockmean; save R, not R_star
                    meanR += (R_blkmean * nInd - meanR) * iIter #blockmean; save R, not R_star 
                    meanR2 += (R2 - meanR2) * iIter
                end
            end

            # #broadcast from rank 0 to other ranks
            # R = MPI.bcast(R, 0, comm)
            # Rinv = inv(R)
        end


        # check convergence
        if iter > burnin

            if sample_delta
                for trait = 1:nTraits
                mcmc_Delta[trait][:, iter - burnin] = deltaArray[trait]
                end
            else
                if iter % outFreq == 0
                    for trait = 1:nTraits
                        mcmc_Delta[trait][:, iout] = deltaArray[trait]
                    end
                end
            end

            if iter % outFreq == 0
                # marker effects
                # for trait = 1:nTraits
                #     mcmc_Alpha[trait][:, iout] = alphaArray[trait]
                # end

                if my_rank == 0
                    # println("iter $iter")
                    # println("nQTL: ", nQTL_sum)
                    # println("Pi count:")
                    # nLoci_array_vec_sum = unflatten(nLoci_array_vec_flatten_sum, length(startPi))
                    # for cat = 1:nCategory
                    #     println("Category $cat")
                    #     nLoci_array_vec_sum_c = nLoci_array_vec_sum[cat]
                    #     iCategori = 1
                    #     for i in keys(startPi)
                    #         println("Pi[$i]: ", nLoci_array_vec_sum_c[iCategori])
                    #         iCategori = iCategori + 1
                    #     end  
                    # end               
                    for cat = 1:nCategory
                        writedlm(outfile["pi"], Pi[cat], ',')
                        writedlm(outfile["beta_effects_variance"], A_vec[cat], ',')
                        writedlm(outfile["genetic_effects_variance"], G_vec_sum[cat], ',')
                        writedlm(outfile["marker_effects_variance"], Atrue_vec_sum[cat] / nQTL_sum[cat], ',')
                    end

                end
                iout += 1
            end
        end

    end # end MCMC iteration loop

    # close outfile
    if my_rank == 0
        for i in keys(outfile)
            close(outfile[i])
        end
    end


    # save mcmc results & posterior mean
    if my_rank == 0  
        for cat in 1:nCategory
            # marker effect variance 
            writedlm(analysis_path * "estA" * string(cat) * ".txt", meanA[cat])
            writedlm(analysis_path * "estA_std" * string(cat) * ".txt", sqrt.((meanA2[cat] .- (meanA[cat] .^ 2))))
            meanAcor[cat] = mean(mcmcAtruecor_c[:, cat][.!isnan.(mcmcAtruecor_c[:, cat])])
            meanAcor2[cat] = mean(mcmcAtruecor_c[:, cat][.!isnan.(mcmcAtruecor_c[:, cat])].^2)
            # beta effect variance
            writedlm(analysis_path * "estB" * string(cat) * ".txt", meanB[cat])
            meanBcor[cat] = mean(mcmcBcor_c[:, cat][.!isnan.(mcmcBcor_c[:, cat])])
            if estimate_vara == true 
                writedlm(analysis_path * "estB_std" * string(cat) * ".txt", sqrt.((meanB2[cat] .- (meanB[cat] .^ 2))))
                meanBcor2[cat] = mean(mcmcBcor_c[:, cat][.!isnan.(mcmcBcor_c[:, cat])] .^ 2)
            end
        end
        writedlm(analysis_path * "estAcor.txt", meanAcor)
        writedlm(analysis_path * "estAcor_std.txt", sqrt.(meanAcor2 .- (meanAcor .^ 2)))
        writedlm(analysis_path * "estBcor.txt", meanBcor)
        if estimate_vara == true
            writedlm(analysis_path * "estBcor_std.txt", sqrt.(meanBcor2 .- (meanBcor .^ 2)))
        end
        
        # genetic variance components
        for cat in 1:nCategory
            writedlm(analysis_path * "estG" * string(cat) * ".txt", meanG[cat])
            writedlm(analysis_path * "estG_std" * string(cat) * ".txt", sqrt.(meanG2[cat] .- (meanG[cat] .^ 2)))
            meanGcor[cat] = mean(mcmcGcor_c[:, cat][.!isnan.(mcmcGcor_c[:, cat])])
            meanGcor2[cat] = mean(mcmcGcor_c[:, cat][.!isnan.(mcmcGcor_c[:, cat])] .^ 2)
            
            meanGenr_ssq[cat] = mean(mcmcGenr_ssq_c[:, cat][.!isnan.(mcmcGenr_ssq_c[:, cat])])
            meanGenr2_ssq[cat] = mean(mcmcGenr_ssq_c[:, cat][.!isnan.(mcmcGenr_ssq_c[:, cat])] .^ 2)
            meanGenr_whq[cat] = mean(mcmcGenr_whq_c[:, cat][.!isnan.(mcmcGenr_whq_c[:, cat])])
            meanGenr2_whq[cat] = mean(mcmcGenr_whq_c[:, cat][.!isnan.(mcmcGenr_whq_c[:, cat])] .^ 2)
        end
        writedlm(analysis_path * "estGcor.txt", meanGcor)
        writedlm(analysis_path * "estGcor_std.txt", sqrt.(meanGcor2 .- (meanGcor .^ 2)))
        writedlm(analysis_path * "estGenr_ssq.txt", meanGenr_ssq)
        writedlm(analysis_path * "estGenr_ssq_std.txt", sqrt.(meanGenr2_ssq .- (meanGenr_ssq .^ 2)))
        writedlm(analysis_path * "estGenr_whq.txt", meanGenr_whq)
        writedlm(analysis_path * "estGenr_whq_std.txt", sqrt.(meanGenr2_whq .- (meanGenr_whq .^ 2)))
        
        #  total genetic variance
        writedlm(analysis_path * "estGtotal.txt", meanGtotal, ',')
        writedlm(analysis_path * "estGtotal_std.txt", sqrt.((meanGtotal2 .- (meanGtotal .^ 2))), ',')
        meanGcor_total = mean(mcmcGcor_total[.!isnan.(mcmcGcor_total)])
        meanGcor_total2 = mean(mcmcGcor_total[.!isnan.(mcmcGcor_total)] .^ 2)
        writedlm(analysis_path * "estGcor_total.txt", meanGcor_total)
        writedlm(analysis_path * "estGcor_total_std.txt", sqrt(meanGcor_total2 - (meanGcor_total^2)))

        # save mcmcGcov_c and mcmcGcov_total
        writedlm(analysis_path * "mcmcGcov_c.txt", mcmcGcov_c)
        writedlm(analysis_path * "mcmcGcov_total.txt", mcmcGcov_total)

        # save mcmcGcor_c and mcmcGcor_total
        writedlm(analysis_path * "mcmcGcor_c.txt", mcmcGcor_c)
        writedlm(analysis_path * "mcmcGcor_total.txt", mcmcGcor_total)

        # save mcmcAtruecor_c
        writedlm(analysis_path * "mcmcAtruecor_c.txt", mcmcAtruecor_c)

        if estimate_vare == true
            writedlm(analysis_path * "estR.txt", meanR)
            writedlm(analysis_path * "estR_std.txt", sqrt.((meanR2 .- (meanR .^ 2))))
        end
        if estimate_pi == true
            for cat in 1:nCategory
                writedlm(analysis_path * "estPi" * string(cat) * ".txt", mean_pi[cat])
                std_pi = deepcopy(mean_pi[cat])
                for i in keys(mean_pi[cat])
                    std_pi[i] = sqrt(mean_pi2[cat][i] - mean_pi[cat][i]^2)
                end
                writedlm(analysis_path * "estPi_std" * string(cat) * ".txt", std_pi)
            end
        end
    end

    writedlm(analysis_path * "mcmc_Delta1.rank$my_rank.txt", mcmc_Delta[1])
    writedlm(analysis_path * "mcmc_Delta2.rank$my_rank.txt", mcmc_Delta[2])
    writedlm(analysis_path * "meanAlpha1.rank$my_rank.txt", meanAlpha[1])
    writedlm(analysis_path * "meanAlpha2.rank$my_rank.txt", meanAlpha[2])

    if my_rank == 0
        time_end = now()
        time_diff = (time_end - time_start).value / 60000 #milliseconds to min
        println("End time: ", time_end)
        println("Running Time (min): ", time_diff)
    end
end

MPI.Init()
@time runMPI()
MPI.Finalize()

