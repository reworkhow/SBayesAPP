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

        "--annotation_size"
        help = ""
        arg_type = Float64
        default = 0.1

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

        "--trait"
        help = ""
        arg_type = Int64
        default = 1
    end

    return parse_args(s)
end

# Use the parsed arguments
args = parse_commandline()

##############################################################
## Input data
##############################################################
data_folder = "/common/zhao/tianjing/simu_chr1_output_new/"
## Input data
##############################################################

nInd = args["sample_size"]
annotation_size = args["annotation_size"]
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
trait = args["trait"]

data_folder_name = "h2_trait1.$h21.h2_trait2.$(h22)_pleioPercent$(pleio_percent)_sampleSize$(nInd)_annotationSize$(annotation_size)_seed$seed/"
data_path = data_folder * data_folder_name
analysis_path = analysis_folder * data_folder_name * "Trait$trait/"
mkpath(analysis_path)
cd(analysis_path)

function sample_variance(x, n, df, scale)
    return (dot(x, x) + df * scale) / rand(Chisq(n + df))
end

function samplePi(nEffects::Number, nTotal::Number)
    return rand(Beta(nTotal - nEffects + 1, nEffects + 1))
end


############################################################################################################
# Input data 
nBlocks = 49
annotationName = ["category1", "category2"]; # ordered by continuous category then categorical category

SNPc1 = readdlm(data_path * "SNPc1.txt",)[:,1]
nSNPc1 = length(SNPc1)
SNPc2 = readdlm(data_path * "SNPc2.txt",)[:,1]
nSNPc2 = length(SNPc2)
annot_size = [nSNPc1, nSNPc2] # number of SNPs in each category

nCon = 0 # number of continuous annotation
nCat = length(annotationName) # number of categorical annotation
annotationType = repeat(["category"], nCat) # ordered by "continue" then "category"
estimate_vare = true
estimate_vara = true
estimate_pi = true 
############################################################################################################

vare = 1 / nInd # residual variance (assuming vary = 1)

#marker effect variance-covariance matrix
varEffects_vec = fill(0.0, nCat) 
Gprior_vec = [zeros(2, 2) for c in 1:nCat]
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
#     Gprior_vec[c] = Gprior_vec[c] / annot_size[c]

#     varEffects_vec[c] = Gprior_vec[c][trait, trait]
# end

for c in 1:nCat
    marker_var_file_path = data_path * "scaled_pleioQTLeff_cov_c$(c).txt"
    if isfile(marker_var_file_path)
        Gprior_vec[c] = reshape(readdlm(marker_var_file_path), 2, 2)
    else
        error("The file $marker_var_file_path does not exist!!!")
    end

    varEffects_vec[c] = Gprior_vec[c][trait, trait]
end






function runMPI(; 
    startPi=0.99, 
    nIter=nIter, outFreq=outFreq, seed=seed, burnin= Int(0.4 * nIter),
    estimate_vare=estimate_vare, estimate_vara=estimate_vara,
    estimate_pi=estimate_pi,
    nInd=nInd, nMarker=nMarker, nBlocks=nBlocks,
    annotationType=annotationType, annotationName=annotationName,
    nCon=nCon, nCat=nCat,
    analysis_path=analysis_path, data_path=data_path,
    varEffects_vec=varEffects_vec, vare=vare,
    annot_size = annot_size)
 
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
    my_anno_matrix_dict = JLD2.load(data_path * "nrank$cluster_size.eigen/bhatXsj/995Eigen/rank$my_rank.anno_matrix_cell_type_human_total_unoverlap_dict.jld2")["my_anno_matrix_dict"]

    # hyper-parameters for varEffects (marker effect variance) and vare
    nCategory = nCat + nCon
    pi_vec = fill(startPi, nCategory)

    if estimate_vare == true
        nue = 4
        scalee = (nue - 2) / nue * vare
    end

    if my_rank == 0
        if estimate_vara == true
            nub = 4
            scaleb_vec = [(nub - 2) / nub * varEffects_vec[c] for c in 1:nCategory]
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

    betaArray = zeros(my_nsnp * nCategory)
    alphaArray = zeros(my_nsnp * nCategory)
    deltaArray = ones(my_nsnp * nCategory)

    #posterior mean 
    meanAlpha = zeros(my_nsnp * nCategory)


    if my_rank == 0 
        if estimate_pi == true 
            mean_pi = zeros(nCategory)
            mean_pi2 = zeros(nCategory)
        end
        
        if estimate_vara == true 
            # beta effect variance
            mean_vara2 = zeros(nCategory)
        end

              
            
            mean_vara = zeros(nCategory)

            # genetic effect variance
            mean_varg = zeros(nCategory)
            mean_varg2 = zeros(nCategory)

            mean_varg_total = 0.0
            mean_varg_total2 = 0.0

        if estimate_vare == true
            mean_vare = 0.0
            mean_vare2 = 0.0
        end
    end
    
    if my_rank == 0
        # outfile to save mcmc results 
        outfile = Dict{String,IOStream}()
        outvar = ["pi"]
        push!(outvar, "beta_effects_variance")
        push!(outvar, "genetic_effects_variance")
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

    ############################################################################
    # print information
    ############################################################################
    if my_rank == 0
        println("---------------- Summary Start --------------")
        println("nIter=$nIter, outFreq=$outFreq, seed=$seed, burnin = $burnin")
        println("Annotation size: ", annotation_size)
        println("Pleiotropy percent: ", pleio_percent)
        println("Heritability for trait 1 (h21): ", h21)
        println("Heritability for trait 2 (h22): ", h22)
        println("Number of ranks: ", nrank)
        println("estimate_vare=$estimate_vare,estimate_vara=$estimate_vara")
        println("estimate_pi=$estimate_pi")
        println("analysis_path=$analysis_path")
        println("data_path=$data_path")
        println("marker effect variance for C1 is: ", varEffects_vec[1])
        println("marker effect variance for C2 is: ", varEffects_vec[2])
        println("residual variance is: ", round.(vare * nInd, digits=3), "  (R*Ind)")
        println("pi=$startPi")
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
    vare_blk = fill(vare, my_nblk)

    
    @showprogress "running MCMC ..." for iter = 1:nIter
        if iter > burnin
            iIter = 1.0 / (iter - burnin)
        end

        varg_blk_cat = [zeros(nCategory) for _ in 1:my_nblk] # varg for different category (saved for hsq computation)
        # compute total varg for each trait without split into different categories
        totalvarg_blk = fill(0.0, my_nblk) 

        nLoci_vec = fill(0, nCategory) # number of causal SNP in each category in this rank
        # sse_vec Computed as dot(\beta, \beta) for whole rank
        sse_vec = fill(0.0, nCategory)


        #for loop for LD blocks
        for b in 1:my_nblk
            blk = my_blkID[b] #block ID
            xpx_vec = xpx_dict[blk]
            xArray_vec = xArray_dict[blk]
            wArray = my_TransformedY_dict[blk][trait]
            annotationMatb = my_anno_matrix_dict[blk]  #annotation boolean data for current blk, nMarkerb-by-C
            SNPIndexb = my_blkSNPsIndex_dict[blk]
            nEigenb, nMarkerb = size(my_TransformedX_dict[blk]) # q & nsnpb

            # # initialize xArrayc/xpxc 
            xArrayc = xArray_vec[1] # nEigenb x nMarkerb matrix
            xpxc = xpx_vec[1]

            vare = vare_blk[b] # VareDn
            invVareDn = 1.0 / vare
            for cat = 1:nCategory
                annoindexc = annotationMatb[:, cat] # SNPs in/out current category
                if nCon != 0
                    if cat <= nCon + 1 # continuous group + categorical group (after nCon+1, xArrayc/xpxc will not change)
                        xArrayc = xArray_vec[cat]
                        xpxc = xpx_vec[cat]
                    end
                end
                
                varEffects = varEffects_vec[cat]
                invVarEffects = 1 / varEffects
                logVarEffects = log(varEffects)
                
                pi = pi_vec[cat]
                logPi = log(pi)
                logPiComp = log(1 - pi)
                logDelta0 = logPi

                invLhs = 1.0 / (invVareDn + invVarEffects)
                logInvLhsMsigma = log(invLhs) - logVarEffects

                for marker = 1:nMarkerb
                    true_marker_num = SNPIndexb[marker] #marker position across my_nsnp
                    if annoindexc[marker] # skip sampling if the marker is not in the category
                        markerIndex = (cat - 1) * my_nsnp + true_marker_num # exact position in betaArray; betaArray: nMarker*nCategory-by-1
                        x = xArrayc[:, marker]
                        oldAlpha = alphaArray[markerIndex]
                        rhs = (dot(x, wArray) + oldAlpha) * invVareDn
                        uhat = invLhs * rhs
                        #sample δj
                        logDelta1 = 0.5 * (logInvLhsMsigma + uhat * rhs) + logPiComp
                        probDelta1 = 1.0 / (1.0 + exp(logDelta0 - logDelta1))

                        #sample marker effects
                        if (rand() < probDelta1)
                            deltaArray[markerIndex] = 1
                            betaArray[markerIndex] = alphaArray[markerIndex] =  uhat + randn() * sqrt(invLhs)
                            BLAS.axpy!(oldAlpha - alphaArray[markerIndex], x, wArray)
                            #update nloci 
                            nLoci_vec[cat] += 1
                        else
                            betaArray[markerIndex] = randn() * sqrt(varEffects)
                            alphaArray[markerIndex] = 0
                            deltaArray[markerIndex] = 0

                            if (oldAlpha != 0.0)
                                BLAS.axpy!(oldAlpha, x, wArray)
                            end
                    
                        end

                        if iter > burnin
                            meanAlpha[markerIndex] += (alphaArray[markerIndex] - meanAlpha[markerIndex]) * iIter
                        end
                    end # end if loop  
                end # end marker loop
            end  # end annotation loop 

            my_TransformedY_dict[blk][trait] = wArray

            # compute genetic variance and heritability 
            XAb = xArray_vec[1]
            what_array_total = zeros(nEigenb)
            
            for cat in 1:nCategory
                # if nCon != 0
                #     if cat <= nCon + 1 # continuous group + categorical group (after nCon+1, xArrayc/xpxc will not change)
                #         XAb = xArray_vec[cat]
                #     end
                # end
                what_array_c = zeros(nEigenb)
                alphaArray_c = alphaArray[(cat - 1) * my_nsnp .+ SNPIndexb]
                what_array_c[:] = XAb * alphaArray_c
                varg_blk_cat[b][cat] = dot(what_array_c, what_array_c) # genetic variance for different category in bth block

                what_array_total += what_array_c
            end
            # total genetic variance
            totalvarg_blk[b] = dot(what_array_total, what_array_total)
            
            if estimate_vare == true
                #sample vare
                vare_blk[b] = sample_variance(wArray, nEigenb, nue, scalee)
            end
        end # end block loop 
        # summing residual variance
        if estimate_vare == true
            vare_blk_sum_rank = sum(vare_blk)
        end


        ########################
        ### Step2. sample Pi ###
        ########################
        # gather nLoci_vec from all ranks 
        nLoci_vec_sum = MPI.Reduce(nLoci_vec, +, 0, comm) # sum the results from all ranks
        if estimate_pi == true
            if my_rank == 0
                for cat = 1:nCategory
                    pi_vec[cat] = samplePi(nLoci_vec_sum[cat], annot_size[cat])
                    if iter > burnin
                        mean_pi .+= (pi_vec .- mean_pi) * iIter
                        mean_pi2 .+= (pi_vec .^ 2 .- mean_pi2) * iIter
                    end
                end
            end
        end

        # broadcast pi_vec from 0 to other ranks 
        pi_vec = MPI.bcast(pi_vec, 0, comm)
        MPI.Barrier(comm)


        
        ### Step3. sample marker variance ###
       for cat in 1:nCategory
            beta_i = betaArray[((cat-1)*my_nsnp+1):((cat-1)*my_nsnp+my_nsnp)]
            sse_vec[cat] = dot(beta_i, beta_i)
       end
        MPI.Barrier(comm)
        # gather the results from all ranks
        sse_vec_sum = MPI.Reduce(sse_vec, +, 0, comm)
        if my_rank == 0
            for cat = 1:nCategory
                if estimate_vara == true
                    varEffects_vec[cat] = (sse_vec_sum[cat] + nub * scaleb_vec[cat]) / rand(Chisq(annot_size[cat] + nub))
                end
                if iter > burnin
                    mean_vara[cat] += (varEffects_vec[cat] - mean_vara[cat]) * iIter
                    if estimate_vara == true
                        mean_vara2[cat] += (varEffects_vec[cat] ^ 2 - mean_vara2[cat]) * iIter
                    end
                    
                end
            end
        end
        if estimate_vara == true
            #broadcast varEffects_vec from 0 to other ranks
            varEffects_vec = MPI.bcast(varEffects_vec, 0, comm)
            MPI.Barrier(comm)
        end

        # summing genetic variance
        if iter > burnin
            varg_cat = sum(varg_blk_cat) # sum varg_blk_cat to get varg for different category across all blocks (vector)
            totalvarg = sum(totalvarg_blk) # sum totalvarg_blk to get total genetic variance across all blocks (value)

            varg_cat_sum = MPI.Reduce(varg_cat, +, 0, comm) # sum the results from all ranks
            totalvarg_sum = MPI.Reduce(totalvarg, +, 0, comm) # sum the results from all ranks

            MPI.Barrier(comm)
            if my_rank == 0
                mean_varg .+= (varg_cat_sum .- mean_varg) * iIter
                mean_varg2 .+= (varg_cat_sum .^ 2 .- mean_varg2) * iIter
                mean_varg_total += (totalvarg_sum - mean_varg_total) * iIter
                mean_varg_total2 += (totalvarg_sum ^ 2 - mean_varg_total2) * iIter
            end
        end
        MPI.Barrier(comm)

        if estimate_vare == true
            # sampling residual variance
            vare_blk_sum = MPI.Reduce(vare_blk_sum_rank, +, 0, comm)
            MPI.Barrier(comm)
            if my_rank == 0
                vare_blkmean = vare_blk_sum / nBlocks
                if iter > burnin
                    vare2 = (vare_blkmean * nInd) ^2
                    mean_vare += (vare_blkmean * nInd - mean_vare) * iIter
                    mean_vare2 += (vare2 - mean_vare2) * iIter
                end
            end
        end

        # check convergence
        if iter > burnin
            if iter % outFreq == 0
                # marker effects
                # for trait = 1:nTraits
                #     mcmc_Alpha[trait][:, iout] = alphaArray[trait]
                # end

                if my_rank == 0
                    println("iter $iter")
                    println("nQTL: ", nLoci_vec_sum)
                    println("pi_vec:", pi_vec)
                    
                    writedlm(outfile["pi"], pi_vec, ',')
                    writedlm(outfile["beta_effects_variance"], varEffects_vec, ',')
                    writedlm(outfile["genetic_effects_variance"], varg_cat_sum, ',')

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
        # beta effect variance
        writedlm(analysis_path * "mean_vara.txt", mean_vara)
        if estimate_vara == true
            writedlm(analysis_path * "std_vara.txt", sqrt.(mean_vara2 .- mean_vara .^ 2))
        end
        # genetic variance components
        writedlm(analysis_path * "mean_varg.txt", mean_varg)
         writedlm(analysis_path * "std_varg.txt", sqrt.(mean_varg2 .- mean_varg .^ 2))
        writedlm(analysis_path * "mean_varg_total.txt", mean_varg_total)
        writedlm(analysis_path * "std_varg_total.txt", sqrt(mean_varg_total2 - mean_varg_total ^ 2))

        if estimate_vare == true
            writedlm(analysis_path * "mean_vare.txt", mean_vare)
            writedlm(analysis_path * "std_vare.txt", sqrt(mean_vare2 - mean_vare ^ 2))            
        end
        if estimate_pi == true
            writedlm(analysis_path * "mean_pi.txt", mean_pi)
            writedlm(analysis_path * "std_pi.txt", sqrt.(mean_pi2 .- mean_pi .^ 2))
        end
    end
        
    writedlm(analysis_path * "meanAlpha.rank$my_rank.txt", meanAlpha)
        

    ############################################################################
    # Step5. show running time
    ############################################################################
    if my_rank == 0
        time_end = now()
        time_diff = (time_end - time_start).value / 60000 #milliseconds to min
        println("End time: ", time_end)
        println("Running Time (min): ", time_diff)
    end
end #end of function

#run MPI
MPI.Init()
@time runMPI()
MPI.Finalize()
