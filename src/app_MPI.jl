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

include(joinpath(@__DIR__, "io", "annotations.jl"))
include(joinpath(@__DIR__, "model", "continuation_state.jl"))
include(joinpath(@__DIR__, "model", "mcmc_setup.jl"))
include(joinpath(@__DIR__, "model", "priors.jl"))
include(joinpath(@__DIR__, "model", "r_blk_state.jl"))
include(joinpath(@__DIR__, "model", "block_setup.jl"))
include(joinpath(@__DIR__, "model", "utilities.jl"))

data_path = ARGS[1] #/common/zhao/jyqqu/MTSBayesCC/data/real_data/T2D_FG/
analysis_path = ARGS[2] #/common/zhao/jyqqu/MTSBayesCC/analysis/real_data/cell_type_human_total_unoverlap/MTSBayesCC_20K/T2D_FG_seed3/
nIter = ARGS[3] # 20K
nIter = parse(Int64, nIter)
seed = ARGS[4] # 3
seed = parse(Int64, seed)
nrank = ARGS[5] # 20
nrank = parse(Int64, nrank)
annot_file = ARGS[6] # ../cell_type_annot_human_total_unoverlap.txt
annot_dict = ARGS[7] # anno_matrix_cell_type_human_total_unoverlap_dict
outFreq = ARGS[8] # 800
outFreq = parse(Int64, outFreq)
starting_value_dir = ARGS[9] # /common/zhao/jyqqu/MTSBayesCC/analysis/real_data/cell_type_human_total_unoverlap/XX/T2D_FG_seed2/
                             # if !is_continue, starting_value_dir is ST_path
                             # folders that save fixed_hyperparameters 
secondary_starting_value_dir = ARGS[10] # folder that contains Gscale files or effects files for split chr scenarios 
                                        # /common/zhao/jyqqu/MTSBayesCC/analysis/real_data/cell_type_human_total_unoverlap/XX/T2D_FG_seed2/
ST_path = ARGS[11]
thin = ARGS[12] # 200
thin = parse(Int64, thin)

N1 = ARGS[13] # average N for Trait 1
N1 = parse(Int64, N1)

N2 = ARGS[14] # average N for Trait 2
N2 = parse(Int64, N2)

has_nCon_input = length(ARGS) ≥ 16 && !(lowercase(ARGS[15]) in ("true", "false"))
nCon_input = has_nCon_input ? parse(Int64, ARGS[15]) : 0
estimate_pi_index = has_nCon_input ? 16 : 15

estimate_pi = ARGS[estimate_pi_index] # estimate pi
# convert estimate_pi to boolean
estimate_pi = (estimate_pi == "true") ? true : false

fixed_hyperparameters = length(ARGS) ≥ estimate_pi_index + 1 ? ARGS[estimate_pi_index + 1] : "false"
fixed_hyperparameters = (fixed_hyperparameters == "true") ? true : false

is_continue = length(ARGS) ≥ estimate_pi_index + 2 ? ARGS[estimate_pi_index + 2] : "true"
is_continue = (is_continue == "true") ? true : false

chr = length(ARGS) ≥ estimate_pi_index + 3 ? ARGS[estimate_pi_index + 3] : ""

if !isempty(chr)
    Split_chr = true
    println("Analysis is done by splitting chromosomes, Chr = $chr")
    data_path = data_path * "SplitByChr/Chr$chr/"
    annot_file = "../../" * annot_file
    analysis_path = analysis_path * "Chr$chr/"
else
    Split_chr = false
    println("Analysis is done by using whole genome data")
end

mkpath(analysis_path)
cd(analysis_path)

# print out N1 and N2
println("average N1: $N1")
println("average N2: $N2")

function sample_variance_sumstats(ycorr_array, nobs, df, scale)
    ntraits = length(ycorr_array)
    SSE = zeros(ntraits, ntraits)
    for traiti = 1:ntraits
        ycorri = ycorr_array[traiti]
        for traitj = traiti:ntraits
            ycorrj = ycorr_array[traitj]
            # multiple the nInd for each trait
            SSE[traiti, traitj] =  dot(ycorri, ycorrj) 
            SSE[traitj, traiti] = SSE[traiti, traitj]
        end
    end
    R = rand(InverseWishart(df + nobs, convert(Array, Symmetric(scale + SSE))))
    return R
end

function invgamma_rand(df, scale)
    return 1 / rand(Gamma(df, 1 / scale))
end

function sample_variance_diag(ycorr_array, nobs, df, scale, nind_array)
    ntraits = length(ycorr_array)
    sampled_diag = zeros(ntraits)
    for traiti = 1:ntraits
        ycorri = ycorr_array[traiti]
        SSE = dot(ycorri, ycorri) * nind_array[traiti]
        scalei = scale[traiti, traiti]
        # This samples from inverse-gamma ~ InvChiSq(df + nobs, scale + SSE)
        sampled_diag[traiti] = invgamma_rand((df + nobs) / 2, (scalei + SSE) / 2)    
    end
    return Diagonal(sampled_diag)
end

############################################################################################################
# Input data 
nMarker = 1154522
nBlocks = 591
# Move this part to function my_rank == 0 
annotation_metadata = load_annotation_metadata(data_path, annot_file; nCon=nCon_input)
annotationName = annotation_metadata.annotationName # ordered by continuous category then categorical category
nLoci_annot = annotation_metadata.nLoci_annot
nCon = annotation_metadata.nCon # number of continuous annotation
nCat = annotation_metadata.nCat # number of categorical annotation
annotationType = annotation_metadata.annotationType # ordered by "continue" then "category"
estimate_pi = estimate_pi
estimate_vare = true
estimate_vara = true
estimate_Gscale = true
estGscale_iter = 2000

ST_version = false
if ST_version
    println("ST version is used.")
end 

if fixed_hyperparameters
    # all the hyper-parameters are fixed in chromosome splitting analysis
    estimate_vare = false
    estimate_vara = false
    estimate_pi = false
    println("All hyper-parameters are fixed.")
else
    if Split_chr
        estimate_vare = false
        println("Residual variance is fixed in chromosome splitting analysis.")
    end
end
############################################################################################################
start_pi_result = build_start_pi(ST_path; estimate_pi=estimate_pi)
startPi = start_pi_result.startPi
Pi00 = start_pi_result.Pi00
if !estimate_pi
    println("Pi is fixed as $startPi and not estimated in the analysis.")
end

Gprior_vec = build_gprior_vec(ST_path, nLoci_annot, Pi00)

function runMPI(; 
    startPi=startPi,
    nIter=nIter, outFreq=outFreq, seed=seed,
    estimate_vare=estimate_vare, estimate_vara=estimate_vara,
    estimate_pi=estimate_pi,
    nMarker=nMarker, nBlocks=nBlocks,
    annotationType=annotationType, annotationName=annotationName,
    nCon=nCon, nCat=nCat,
    analysis_path=analysis_path, data_path=data_path,
    Gprior_vec=Gprior_vec,
    thin=thin, N1 = N1, N2 = N2,
    estimate_Gscale=estimate_Gscale, estGscale_iter=estGscale_iter,
    save_delta = true)

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
    my_nGWAS_dict = JLD2.load(data_path * "nrank$cluster_size.eigen/bhatXsj/995Eigen/rank$my_rank.nGWAS_dict.jld2")["my_nGWAS_dict"]
    my_nblk = length(my_blkID)
    my_nsnp = sum(map(length, values(my_blkSNPsIndex_dict)))
    # need to change accordingly
    my_anno_matrix_dict = JLD2.load(data_path * "nrank$cluster_size.eigen/bhatXsj/995Eigen/rank$my_rank.$annot_dict.jld2")["my_anno_matrix_dict"] 

    if is_continue
        burnin = 0
    else 
        burnin = 1000 
    end

    # hyper-parameters for A (marker effect variance) and R
    nCategory = nCat + nCon
    # marker effect variance matrix
    A_vec = [zeros(2, 2) for c in 1:nCategory]
    if fixed_hyperparameters
        A_vec_starting_path = starting_value_dir
        for c in 1:nCategory
            A_vec[c] = readdlm(A_vec_starting_path * "estB$c.txt")
        end
    else
        if is_continue
            A_vec_starting_path = starting_value_dir * "beta_effect_var_matrices_last_sample/"
            for c in 1:nCategory
                A_vec[c] = readdlm(A_vec_starting_path * "beta_effect_matrix_$c.txt", ',')
            end
        else # starting values 
            A_vec = Gprior_vec
        end
    end
    Ainv_vec = [inv(A_vec[c]) for c in 1:nCategory]

    if fixed_hyperparameters
        Pi_starting_path = starting_value_dir
        Pi = [Dict{Vector{Float64},Float64}() for c in 1:nCategory]
        for c in 1:nCategory
            Pi[c] = read_to_dict_posterior_mean(Pi_starting_path * "estPi$c.txt")
        end
    else
        if is_continue
            Pi_starting_path = starting_value_dir * "pi_last_sample/"
            Pi = [Dict{Vector{Float64},Float64}() for c in 1:nCategory]
            for c in 1:nCategory
                Pi[c] = read_to_dict(Pi_starting_path * "pi_$c.txt")
            end
        else
            Pi = [deepcopy(startPi) for c in 1:nCategory]
        end
    end

    Rprior = [1. / N1 0. ; 0. 1. / N2]
    nTraits = 2

    if estimate_vare == true
        df_R = 4 + nTraits
        scale_R = Rprior * (df_R - nTraits - 1)
    end

    if my_rank == 0
        if estimate_vara == true
            df_G = 4 + nTraits

            if is_continue && estimate_Gscale
                scale_G_vec = [zeros(2, 2) for c in 1:nCategory]
                for c in 1:nCategory
                    scale_G_vec[c] = readdlm(secondary_starting_value_dir * "scale_G$c.txt")
                end
                estimate_Gscale = false # estimate_Gscale = false for continuous analysis
                println("scale_G_vec is fixed as saved in $secondary_starting_value_dir")
            else
                scale_G_vec = [Gprior_vec[c] * (df_G - nTraits - 1) for c in 1:nCategory]
                println("scale_G_vec is computed by ST h2.")
            end
        end
    end

    # Broadcast `estimate_Gscale` from rank 0 to all ranks
    MPI.Barrier(comm)
    if is_continue && estimate_Gscale
        estimate_Gscale = false
        println("estimate_Gscale is set as $estimate_Gscale for all ranks!")
    end
    MPI.Barrier(comm)

    if my_rank == 0
        writedlm(analysis_path * "annotationName.txt", annotationName)
    end
    xpx_dict, xArray_dict, my_anno_matrix_dict = prepare_block_state!(
        my_blkID,
        my_blkSNPsIndex_dict,
        my_TransformedX_dict,
        my_anno_matrix_dict,
        nCon,
        annotationType,
    )

    #output
    β = zeros(nTraits)
    newα = zeros(nTraits)  #α=Dβ
    oldα = zeros(nTraits)
    w = zeros(nTraits)
    δ = zeros(nTraits)

    if !is_continue
        betaArray = [zeros(my_nsnp * nCategory) for t in 1:nTraits] #-> ordered by annotation groups 
        alphaArray = [zeros(my_nsnp * nCategory) for t in 1:nTraits]
        deltaArray = [zeros(my_nsnp * nCategory) for t in 1:nTraits]
    else
        if Split_chr
            effect_starting_path = secondary_starting_value_dir * "Chr$chr/"
            delta_starting_path = secondary_starting_value_dir * "Chr$chr/last_sample_delta/"
        else
            effect_starting_path = starting_value_dir
            delta_starting_path = starting_value_dir * "last_sample_delta/"
        end
        betaArray, alphaArray, deltaArray = initialize_effect_state!(
            effect_starting_path,
            delta_starting_path,
            my_rank,
            my_nsnp,
            nCategory,
            nTraits,
            my_blkID,
            my_blkSNPsIndex_dict,
            my_TransformedX_dict,
            my_TransformedY_dict,
        )
    end

    #posterior mean 
    meanAlpha = [zeros(my_nsnp * nCategory) for t in 1:nTraits]
    
    nOutput = Int(floor((nIter-burnin) / outFreq))
    if save_delta
        mcmc_Delta = [zeros(my_nsnp * nCategory, nOutput) for t in 1:nTraits]
    end

    if my_rank == 0
        rank0_mcmc_state = initialize_rank0_mcmc_state(
            Pi,
            nIter,
            burnin,
            thin,
            nTraits,
            nCategory;
            estimate_pi=estimate_pi,
            estimate_vara=estimate_vara,
            estimate_vare=estimate_vare,
        )
        nsample4mean = rank0_mcmc_state.nsample4mean
        mean_pi = rank0_mcmc_state.mean_pi
        mean_pi2 = rank0_mcmc_state.mean_pi2
        meanB2 = rank0_mcmc_state.meanB2
        meanA2 = rank0_mcmc_state.meanA2
        meanBcor2 = rank0_mcmc_state.meanBcor2
        meanAcor2 = rank0_mcmc_state.meanAcor2
        meanA = rank0_mcmc_state.meanA
        meanAcor = rank0_mcmc_state.meanAcor
        meanB = rank0_mcmc_state.meanB
        meanBcor = rank0_mcmc_state.meanBcor
        meanG = rank0_mcmc_state.meanG
        meanG2 = rank0_mcmc_state.meanG2
        meanGcor = rank0_mcmc_state.meanGcor
        meanGcor2 = rank0_mcmc_state.meanGcor2
        meanSSE = rank0_mcmc_state.meanSSE
        meanGtotal = rank0_mcmc_state.meanGtotal
        meanGtotal2 = rank0_mcmc_state.meanGtotal2
        mcmcAtruecor_c = rank0_mcmc_state.mcmcAtruecor_c
        mcmcBcor_c = rank0_mcmc_state.mcmcBcor_c
        mcmcGcov_c = rank0_mcmc_state.mcmcGcov_c
        mcmcGcor_c = rank0_mcmc_state.mcmcGcor_c
        mcmcGcov_total = rank0_mcmc_state.mcmcGcov_total
        mcmcGcor_total = rank0_mcmc_state.mcmcGcor_total
        meanR = rank0_mcmc_state.meanR
        meanR2 = rank0_mcmc_state.meanR2
    end

    file_names = nothing
    if my_rank == 0
        file_names = prepare_mcmc_output_files(analysis_path)
    end

    if my_rank == 0
        println("---------------- Summary Start --------------")
        println("nIter=$nIter, outFreq=$outFreq, seed=$seed, burnin = $burnin")
        println("startPi is: $Pi")
        println("Number of ranks: ", nrank)
        println("estimate_vare=$estimate_vare,estimate_vara=$estimate_vara")
        println("estimate_pi=$estimate_pi")
        println("estimate_Gscale=$estimate_Gscale")
        println("analysis_path=$analysis_path")
        println("data_path=$data_path")
        if estimate_vara
            println("scale_G_vec is: ", scale_G_vec)
        end
        println("nMarker=$nMarker, nBlocks=$nBlocks")
        println("nCat = $nCat, nCon = $nCon")
        println("thin = $thin")
        time_start = now()
        println("Start time: ", time_start)
        println("---------------- Summary End ----------------")
    end
    MPI.Barrier(comm)
    println("In rank$my_rank, there are $my_nblk LD blocks, and $my_nsnp SNPs in total.")
    MPI.Barrier(comm)

    nlabel = 4 

    iout = 1
    iter_after_burnin_thin_index = 1
    iter_after_Greset_thin_index = 0
    last_saved_iter = nIter - (nIter - burnin) % outFreq

    R_blk = initialize_r_blk_state(
        my_nblk,
        my_rank,
        Rprior;
        is_continue=is_continue,
        starting_value_dir=starting_value_dir,
        fixed_r_dir=estimate_vare ? nothing : secondary_starting_value_dir,
    )

    estGscale_iter = min(estGscale_iter, nIter)

    # varg Computed as dot(X\alpha, X\alpha)  
    varg_blk_cat = [zeros(nTraits, nCategory) for _ in 1:my_nblk] # varg for different category & trait (saved for hsq computation)
    varg_cov_blk_cat = [zeros(nCategory) for _ in 1:my_nblk] # genetic covariance for different category
    # compute total varg for each trait without split into different categories
    totalvarg_blk = [zeros(nTraits, nTraits) for _ in 1:my_nblk]

    # ssq_blk_cat Computed as dot(\alpha, \alpha)
    ssq_blk_cat = [zeros(nTraits, nCategory) for _ in 1:my_nblk] # sum of square for different category & trait (saved for enrichment)
    ssq_cov_blk_cat = [zeros(nCategory) for _ in 1:my_nblk] # sum of square for different category (saved for enrichment)

    nLoci_array_vec = [fill(0, nlabel) for c in 1:nCategory] # number of loci (pi) in each category for SNPs in this rank
    # SSE_vec Computed as dot(\beta, \beta) for whole rank
    SSE_vec = [zeros(nTraits, nTraits) for c in 1:nCategory]
    # SSQ_vec Computed as dot(\alpha, \alpha) for whole rank -> used for sampling Gscale
    # SSQ_vec = [zeros(nTraits, nTraits) for c in 1:nCategory] 
    # G_vec Computed as dot(X\alpha, X\alpha) for whole rank
    G_vec = [zeros(nTraits, nTraits) for c in 1:nCategory] # genetic variance for different category for SNPs in this rank

    R_blkmean = zeros(nTraits, nTraits)

    @showprogress "running MCMC ..." for iter = 1:nIter
        
        within_burnin_estGscale = estimate_Gscale && iter <= estGscale_iter
        do_thin = iter > burnin && (iter - burnin) % thin == 0

        if do_thin
            iIter = 1.0 / iter_after_burnin_thin_index
        end

        if within_burnin_estGscale
            iIter_scaleG = 1.0 / iter
        end

        # Zero-out per-block arrays
        for b in 1:my_nblk
            fill!(varg_blk_cat[b], 0.0)
            fill!(varg_cov_blk_cat[b], 0.0)
            fill!(totalvarg_blk[b], 0.0)
            fill!(ssq_blk_cat[b], 0.0)
            fill!(ssq_cov_blk_cat[b], 0.0)
        end

        # Zero-out per-category arrays
        for c in 1:nCategory
            fill!(nLoci_array_vec[c], 0.0)
            fill!(SSE_vec[c], 0.0)
            #fill!(SSQ_vec[c], 0.0)
            fill!(G_vec[c], 0.0)
        end

        #for loop for LD blocks
        for b in 1:my_nblk
            blk = my_blkID[b] #block ID
            xpx_vec = xpx_dict[blk]
            xArray_vec = xArray_dict[blk]
            wArray = my_TransformedY_dict[blk]
            annotationMatb = my_anno_matrix_dict[blk]  #annotation boolean data for current blk, nMarkerb-by-C
            SNPIndexb = my_blkSNPsIndex_dict[blk]
            nEigenb, nMarkerb = size(my_TransformedX_dict[blk]) # q & nsnpb
            nInd = my_nGWAS_dict[blk] # n

            # # initialize xArrayc/xpxc 
            xArrayc = xArray_vec[1] # nEigenb x nMarkerb matrix
            xpxc = xpx_vec[1]

            Rinv = inv(R_blk[b])

            MarkerOrder = shuffle(1:nMarkerb)
            for marker = MarkerOrder
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

                            d0 = deepcopy(δ)
                            d1 = deepcopy(δ)
                            d0[k] = 0.0
                            d1[k] = 1.0

                            #sample δj
                            logDelta0 = -0.5 * (log(Ginv11) - gHat0^2 * Ginv11) + log(BigPi[d0]) #logPi
                            logDelta1 = -0.5 * (log(C11) - gHat1^2 * C11) + log(BigPi[d1]) #logPiComp
                            probDelta1 = 1.0 / (1.0 + exp(logDelta0 - logDelta1))

                            #sample marker effects
                            if (rand() < probDelta1) #δj=1
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
                            if do_thin
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
                if ST_version
                    sampled_R = sample_variance_diag(wArray, nEigenb, df_R, scale_R, nInd)
                    for traiti = 1:nTraits
                        thres = sum(ssq_blk_cat[b][traiti, :]) / totalvarg_blk[b][traiti, traiti]
                        if (thres > 1.1)
                            R_blk[b][traiti, traiti] = sampled_R[traiti, traiti]
                        else
                            R_blk[b][traiti, traiti] = 1.0  
                        end
                    end
                else
                    sampled_R = sample_variance_sumstats(wArray, nEigenb, df_R, scale_R)
                    Rcor = compute_correlation(sampled_R)
                    for traiti = 1:nTraits
                        thres = sum(ssq_blk_cat[b][traiti, :]) / totalvarg_blk[b][traiti, traiti]
                        if (thres > 1.1)
                            R_blk[b][traiti, traiti] = sampled_R[traiti, traiti]
                        else
                            R_blk[b][traiti, traiti] = 1.0 / nInd[traiti]
                        end
                    end
                    # tune covariance in R_blk
                    Rcov = Rcor * sqrt(R_blk[b][1, 1] * R_blk[b][2, 2])
                    R_blk[b][1, 2] = R_blk[b][2, 1] = Rcov
                end
            end 
        end # end block loop 

        # summing residual variance
        if estimate_vare == true
            R_blk_sum_rank = sum(R_blk)
        end

        # get \alpha'\alpha to compute enrichment 
        # if do_thin
        #     # compute \alpha'\alpha for each category
        #     SSQ_vec = [zeros(nTraits, nTraits) for c in 1:nCategory]
        #     ssq_cat = sum(ssq_blk_cat)
        #     ssq_cov_cat = sum(ssq_cov_blk_cat)

        #     for cat = 1:nCategory
        #         SSQ_vec[cat][1, 1] = ssq_cat[1, cat]
        #         SSQ_vec[cat][2, 2] = ssq_cat[2, cat]
        #         SSQ_vec[cat][1, 2] = SSQ_vec[cat][2, 1] = ssq_cov_cat[cat]
        #     end

        #     SSQ_vec_flatten = flatten_matrices(SSQ_vec)
        #     MPI.Barrier(comm)
        #     SSQ_vec_flatten_sum = MPI.Reduce(SSQ_vec_flatten, +, 0, comm)
        #     MPI.Barrier(comm)

        #     if my_rank == 0
        #         SSQ_vec_sum = unflatten_matrices(SSQ_vec_flatten_sum, nTraits, nTraits)
        #     end
        # end

        # get true A matrix by alphaArray
        # use only pleiotropic markers to compute QTL effect variance matrix
        if do_thin 
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
            MPI.Barrier(comm)
            Atrue_vec_flatten_sum = MPI.Reduce(Atrue_vec_flatten, +, 0, comm)
            MPI.Barrier(comm)
            nQTL_sum = MPI.Reduce(nQTL, +, 0, comm)
            MPI.Barrier(comm)
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
                    if estimate_vara == true
                        meanA2[cat] += (Atrue_cat .^ 2 - meanA2[cat]) * iIter
                    end
                    # correlation values
                    mcmcAtruecor_c[iter_after_burnin_thin_index, cat] = compute_correlation(Atrue_cat)
                    # save the Atrue_cat to file
                    open(file_names["marker_effects_variance"], "a") do io
                        writedlm(io, Atrue_cat, ',')
                    end
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
        MPI.Barrier(comm)
        nLoci_array_vec_flatten_sum = MPI.Reduce(nLoci_array_vec_flatten, +, 0, comm) # sum the results from all ranks
        MPI.Barrier(comm)
        if my_rank == 0
            nLoci_array_vec_sum = unflatten(nLoci_array_vec_flatten_sum, nlabel)
        end
        if estimate_pi == true
            tempPi_vec = [zeros(nlabel) for cat = 1:nCategory]
            if my_rank == 0
                for cat = 1:nCategory
                    tempPi_vec[cat] = rand(Dirichlet(nLoci_array_vec_sum[cat] .+ 1))
                    if do_thin
                        tempPi2 = tempPi_vec[cat] .^ 2
                        for (iCategori, i) in enumerate(keys(Pi[cat]))
                            mean_pi[cat][i] += (tempPi_vec[cat][iCategori] - mean_pi[cat][i]) * iIter
                            mean_pi2[cat][i] += (tempPi2[iCategori] - mean_pi2[cat][i]) * iIter
                        end
                    end   
                end
            end

            ############################################################
            # broadcast tempPi_vec from 0 to other ranks 
            tempPi_vec_flatten = flatten(tempPi_vec)
            MPI.Barrier(comm)
            tempPi_vec_flatten = MPI.bcast(tempPi_vec_flatten, 0, comm)
            MPI.Barrier(comm)
            tempPi_vec = unflatten(tempPi_vec_flatten, nlabel)
            for cat = 1:nCategory
                iCategori = 1
                for i in keys(Pi[cat])
                    Pi[cat][i] = tempPi_vec[cat][iCategori] #annotation specific pi
                    iCategori = iCategori + 1
                end
            end
        end
        MPI.Barrier(comm)

        if do_thin
            varg_cat = sum(varg_blk_cat) # sum varg_blk_cat to get varg for different category across blocks in this rank
            varg_cov_cat = sum(varg_cov_blk_cat) # sum varg_cov_blk_cat to get genetic covariance for different category across blocks in this rank
            totalvarg = sum(totalvarg_blk) # sum totalvarg_blk to get total genetic variance across all blocks in this rank

            for cat = 1:nCategory
                G_vec[cat][1, 1] = varg_cat[1, cat]
                G_vec[cat][2, 2] = varg_cat[2, cat]
                G_vec[cat][1, 2] = G_vec[cat][2, 1] = varg_cov_cat[cat]
            end

            G_vec_flatten = flatten_matrices(G_vec)
            # sum the G_vec_flatten from all ranks to root rank
            MPI.Barrier(comm)
            G_vec_flatten_sum = MPI.Reduce(G_vec_flatten, +, 0, comm)
            MPI.Barrier(comm)
            G_total = MPI.Reduce(totalvarg, +, 0, comm)
            MPI.Barrier(comm)

            if my_rank == 0
                G_vec_sum = unflatten_matrices(G_vec_flatten_sum, nTraits, nTraits)

                mcmcGcor_total[iter_after_burnin_thin_index] = compute_correlation(G_total)
                mcmcGcov_total[iter_after_burnin_thin_index] = G_total[1, 2]
                #SSQ_total = sum(SSQ_vec_sum)

                for cat = 1:nCategory
                    # correlation values
                    mcmcGcor_c[iter_after_burnin_thin_index, cat] = compute_correlation(G_vec_sum[cat])
                    mcmcGcov_c[iter_after_burnin_thin_index, cat] = G_vec_sum[cat][1, 2]
                    # two ways to compute enrichment
                    #mcmcGenr_ssq_c[iter_after_burnin_thin_index, cat] = SSQ_vec_sum[cat][1, 2] / SSQ_total[1, 2]
                    #mcmcGenr_whq_c[iter_after_burnin_thin_index, cat] = G_vec_sum[cat][1, 2] / G_total[1, 2]
                    meanG[cat] += (G_vec_sum[cat] - meanG[cat]) * iIter
                    meanG2[cat] += (G_vec_sum[cat] .^ 2 - meanG2[cat]) * iIter
                end
                meanGtotal += (G_total - meanGtotal) * iIter
                meanGtotal2 += (G_total .^ 2 - meanGtotal2) * iIter
            end
        end

        # get SSE_vec to sample A in root and to estGscale
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

        # gather the results from all ranks
        MPI.Barrier(comm)
        SSE_vec_flatten_sum = MPI.Reduce(SSE_vec_flatten, +, 0, comm)
        MPI.Barrier(comm)
        if my_rank == 0
            SSE_vec_sum = unflatten_matrices(SSE_vec_flatten_sum, nTraits, nTraits)
        end

        # summing genetic variance
        if within_burnin_estGscale 
            if my_rank == 0
                for cat = 1:nCategory
                    meanSSE[cat] += (SSE_vec_sum[cat] - meanSSE[cat]) * iIter_scaleG
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

        if within_burnin_estGscale
            if my_rank == 0
                # meanG
                # mean_pi
                Gprior_vec = deepcopy(meanSSE)
                for cat = 1:nCategory
                    Gprior_vec[cat] = meanSSE[cat] / nLoci_annot[cat]
                    if ST_version
                        # Only take the diagonal part and scale it
                        diag_vals = diag(Gprior_vec[cat]) * (df_G - nTraits - 1)
                        scale_G_vec[cat] = Diagonal(diag_vals)
                    else
                        scale_G_vec[cat] = Gprior_vec[cat] * (df_G - nTraits - 1)
                    end
                end

                if iter == estGscale_iter
                    # save the scale_G_vec
                    for cat = 1:nCategory
                        writedlm(analysis_path * "scale_G" * string(cat) * ".txt", scale_G_vec[cat])
                    end
                end
            end
        end

        ### Step3. sample variance ###
        if my_rank == 0
            for cat = 1:nCategory
                if estimate_vara == true
                    if ST_version
                        df_A = df_G + nLoci_annot[cat]
                        A_diag = zeros(nTraits)
                        for traiti in 1:nTraits
                            A_vec_sampler = scale_G_vec[cat][traiti, traiti] + SSE_vec_sum[cat][traiti, traiti]
                            A_diag[traiti] = invgamma_rand(df_A / 2, A_vec_sampler / 2)
                        end
                        A_vec[cat] = Diagonal(A_diag)
                    else
                        A_vec_sampler = convert(Array, Symmetric(scale_G_vec[cat] + SSE_vec_sum[cat]))
                        if is_positive_definite(A_vec_sampler)
                            A_vec[cat] = rand(InverseWishart(df_G + nLoci_annot[cat], A_vec_sampler))
                        else
                            println("A_vec_sampler in iteration $iter for category $cat is not positive definite: ")
                            println(A_vec_sampler)
                        end
                    end
                end
                if do_thin
                    # save mean for beta effect variance 
                    meanB[cat] += (A_vec[cat] - meanB[cat]) * iIter
                    # correlation values
                    mcmcBcor_c[iter_after_burnin_thin_index, cat] = compute_correlation(A_vec[cat])
                    if estimate_vara == true
                        meanB2[cat] += (A_vec[cat] .^ 2 - meanB2[cat]) * iIter
                    end
                end
            end
        end
        if estimate_vara == true 
            # broadcast A_vec from 0 to other ranks
            A_vec_flatten = flatten_matrices(A_vec)
            MPI.Barrier(comm)
            A_vec_flatten = MPI.bcast(A_vec_flatten, 0, comm)
            MPI.Barrier(comm)
            A_vec = unflatten_matrices(A_vec_flatten, nTraits, nTraits)
            Ainv_vec = [inv(A_vec[cat]) for cat in 1:nCategory]
        end
        MPI.Barrier(comm)

        # sampling residual variance
        if estimate_vare
            R_blk_sum = MPI.Reduce(R_blk_sum_rank, +, 0, comm)
            MPI.Barrier(comm)
            if my_rank == 0
                R_blkmean = R_blk_sum / nBlocks
                if do_thin
                    R2 = (R_blkmean) .^ 2
                    meanR += (R_blkmean - meanR) * iIter
                    meanR2 += (R2 - meanR2) * iIter
                end
            end
            # # broadcast R_blkmean from 0 to other ranks
            # R_blkmean = MPI.bcast(R_blkmean, 0, comm)
            # for b in 1:my_nblk
            #     R_blk[b] = R_blkmean
            # end
        end


        # check convergence
        if iter > burnin 
             
            if save_delta
                if iter - burnin > 0 && (iter - burnin) % outFreq == 0
                    for trait = 1:nTraits
                        mcmc_Delta[trait][:, iout] = deltaArray[trait]
                    end
                end
            end

            if iter - burnin > 0 && (iter - burnin) % outFreq == 0

                if my_rank == 0
                    println("iter $iter")
                    for cat = 1:nCategory
                        open(file_names["pi"], "a") do io
                            writedlm(io, Pi[cat], ',')
                        end
                        open(file_names["beta_effects_variance"], "a") do io
                            writedlm(io, A_vec[cat], ',')
                        end
                        open(file_names["genetic_effects_variance"], "a") do io
                            writedlm(io, G_vec_sum[cat], ',')
                        end
                    end
                end
                iout += 1
            end

            if do_thin
                if my_rank == 0
                    open(file_names["total_genetic_effects_variance"], "a") do io
                        writedlm(io, G_total, ',')
                    end
                end
                iter_after_burnin_thin_index += 1
            end

            if iter == last_saved_iter
                for t in 1:nTraits
                    writedlm(analysis_path * "last_mcmc_betaArray$t.rank$my_rank.txt", betaArray[t])
                end
                # Save R_blk
                mkpath(analysis_path * "last_sample_R_blk/")  # Create folder if not exist
                for b in 1:my_nblk
                    writedlm(analysis_path * "last_sample_R_blk/R_blk_b$(b)_rank$(my_rank).txt", R_blk[b])
                end
                # Save beta_effect variance matrices (A_vec)
                mkpath(analysis_path * "beta_effect_var_matrices_last_sample/")
                for c in 1:nCategory
                     writedlm(analysis_path * "beta_effect_var_matrices_last_sample/beta_effect_matrix_$(c).txt", A_vec[c], ',')
                end
                # Save Pi from last sample
                mkpath(analysis_path * "pi_last_sample/")
                for c in 1:nCategory
                    open(analysis_path * "pi_last_sample/pi_$(c).txt", "w") do io
                        writedlm(io, Pi[c], ',')
                    end
                end
                ### Save delta (last column of mcmc_Delta)
                if save_delta
                    mkpath(analysis_path * "last_sample_delta/")
                    for t in 1:nTraits
                        last_delta = mcmc_Delta[t][:, end]  # extract last column
                        open(analysis_path * "last_sample_delta/last_sample_delta$(t)_rank$(my_rank).txt", "w") do io
                            writedlm(io, last_delta)
                        end
                    end
                end
            end
        end

        # GC 
        if iter % 1000 == 0
            GC.gc()
        end
    end # end MCMC iteration loop
    

    # save mcmc results & posterior mean
    if my_rank == 0
        for cat in 1:nCategory
            # marker effect variance 
            writedlm(analysis_path * "estA" * string(cat) * ".txt", meanA[cat])
            meanAcor[cat] = mean(mcmcAtruecor_c[:, cat][.!isnan.(mcmcAtruecor_c[:, cat])])
            # beta effect variance
            writedlm(analysis_path * "estB" * string(cat) * ".txt", meanB[cat])
            meanBcor[cat] = mean(mcmcBcor_c[:, cat][.!isnan.(mcmcBcor_c[:, cat])])
            if estimate_vara == true 
                writedlm(analysis_path * "estB_std" * string(cat) * ".txt", sqrt.((meanB2[cat] .- (meanB[cat] .^ 2))))
                writedlm(analysis_path * "estA_std" * string(cat) * ".txt", sqrt.((meanA2[cat] .- (meanA[cat] .^ 2))))
                meanBcor2[cat] = mean(mcmcBcor_c[:, cat][.!isnan.(mcmcBcor_c[:, cat])] .^ 2)
                meanAcor2[cat] = mean(mcmcAtruecor_c[:, cat][.!isnan.(mcmcAtruecor_c[:, cat])].^2)
            end
        end

        # save mcmcGcov_c and mcmcGcov_total
        writedlm(analysis_path * "mcmcGcov_c.txt", mcmcGcov_c)
        writedlm(analysis_path * "mcmcGcov_total.txt", mcmcGcov_total)

        # save mcmcGcor_c and mcmcGcor_total
        writedlm(analysis_path * "mcmcGcor_c.txt", mcmcGcor_c)
        writedlm(analysis_path * "mcmcGcor_total.txt", mcmcGcor_total)

        # save mcmcAtruecor_c
        writedlm(analysis_path * "mcmcAtruecor_c.txt", mcmcAtruecor_c)

        writedlm(analysis_path * "estAcor.txt", meanAcor)
        writedlm(analysis_path * "estBcor.txt", meanBcor)
        if estimate_vara == true
            #writedlm(analysis_path * "estBcor_std.txt", sqrt.(meanBcor2 .- (meanBcor .^ 2)))
            #writedlm(analysis_path * "estAcor_std.txt", sqrt.(meanAcor2 .- (meanAcor .^ 2)))
            writedlm(analysis_path * "meanBcor2.txt", meanBcor2)
            writedlm(analysis_path * "meanAcor2.txt", meanAcor2)
        end

        # genetic variance components
        for cat in 1:nCategory
            writedlm(analysis_path * "estG" * string(cat) * ".txt", meanG[cat])
            writedlm(analysis_path * "estG_std" * string(cat) * ".txt", sqrt.(meanG2[cat] .- (meanG[cat] .^ 2)))
            meanGcor[cat] = mean(mcmcGcor_c[:, cat][.!isnan.(mcmcGcor_c[:, cat])])
            meanGcor2[cat] = mean(mcmcGcor_c[:, cat][.!isnan.(mcmcGcor_c[:, cat])] .^ 2)
            
            #meanGenr_ssq[cat] = mean(mcmcGenr_ssq_c[:, cat][.!isnan.(mcmcGenr_ssq_c[:, cat])])
            #meanGenr2_ssq[cat] = mean(mcmcGenr_ssq_c[:, cat][.!isnan.(mcmcGenr_ssq_c[:, cat])] .^ 2)
            #meanGenr_whq[cat] = mean(mcmcGenr_whq_c[:, cat][.!isnan.(mcmcGenr_whq_c[:, cat])])
            #meanGenr2_whq[cat] = mean(mcmcGenr_whq_c[:, cat][.!isnan.(mcmcGenr_whq_c[:, cat])] .^ 2)
        end
        writedlm(analysis_path * "estGcor.txt", meanGcor)
        writedlm(analysis_path * "estGcor_std.txt", sqrt.(meanGcor2 .- (meanGcor .^ 2)))
        #writedlm(analysis_path * "estGenr_ssq.txt", meanGenr_ssq)
        #writedlm(analysis_path * "estGenr_ssq_std.txt", sqrt.(meanGenr2_ssq .- (meanGenr_ssq .^ 2)))
        #writedlm(analysis_path * "estGenr_whq.txt", meanGenr_whq)
        #writedlm(analysis_path * "estGenr_whq_std.txt", sqrt.(meanGenr2_whq .- (meanGenr_whq .^ 2)))
        
        #  total genetic variance
        writedlm(analysis_path * "estGtotal.txt", meanGtotal, ',')
        writedlm(analysis_path * "estGtotal_std.txt", sqrt.((meanGtotal2 .- (meanGtotal .^ 2))), ',')
        meanGcor_total = mean(mcmcGcor_total[.!isnan.(mcmcGcor_total)])
        meanGcor_total2 = mean(mcmcGcor_total[.!isnan.(mcmcGcor_total)] .^ 2)
        writedlm(analysis_path * "estGcor_total.txt", meanGcor_total)
        writedlm(analysis_path * "estGcor_total_std.txt", sqrt(meanGcor_total2 - (meanGcor_total^2)))

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

    if save_delta
        writedlm(analysis_path * "mcmc_Delta1.rank$my_rank.txt", mcmc_Delta[1])
        writedlm(analysis_path * "mcmc_Delta2.rank$my_rank.txt", mcmc_Delta[2])
    end
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

