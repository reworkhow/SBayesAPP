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
trait = args["trait"]

# data_folder_name = "h2_trait1.$h21.h2_trait2.$(h22)_pleioPercent$(pleio_percent)_sampleSize$(nInd)_annotationSize$(annotation_size)_seed$seed/"
data_folder_name = "h2_trait1.$h21.h2_trait2.$(h22)_pleioPercent$(pleio_percent)_sampleSize$(nInd)_seed$seed/"
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
estimate_vare = true
estimate_vara = true
estimate_pi = true
############################################################################################################

vare = 1 / nInd # residual variance (assuming vary = 1)

#marker effect variance-covariance matrix
# Gprior = zeros(2, 2)
# saved_sdbv = readdlm(data_path * "sdbv.txt")[:, 1]
# saved_sdy = readdlm(data_path * "sdy.txt")[:, 1]

# bv_var_file_path = data_path * "bv_cov_var.txt"
# saved_bv_var = readdlm(bv_var_file_path)
# t1_scaler = 1 / saved_sdbv[1] * sqrt(h21) / saved_sdy[1]
# t2_scaler = 1 / saved_sdbv[2] * sqrt(h22) / saved_sdy[2]
# Gprior[1,1] = saved_bv_var[1,1] * t1_scaler^2
# Gprior[2,2] = saved_bv_var[2,2] * t2_scaler^2
# Gprior[1,2] = Gprior[2,1] = saved_bv_var[1,2] * t1_scaler * t2_scaler
# Gprior = saved_bv_var
# Gprior = Gprior / (nMarker * 0.01)

startPi = 0.9
varEffects = 0.1 / (nMarker * (1-startPi)) # marker effect variance






function runMPI(; 
    startPi=startPi,
    nIter=nIter, outFreq=outFreq, seed=seed, burnin= Int(0.4 * nIter),
    estimate_vare=estimate_vare, estimate_vara=estimate_vara,
    estimate_pi=estimate_pi,
    nInd=nInd, nMarker=nMarker, nBlocks=nBlocks,
    analysis_path=analysis_path, data_path=data_path,
    varEffects=varEffects, vare=vare)
 
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
    
    # hyper-parameters for varEffects (marker effect variance) and vare
    
    pi = startPi

    if estimate_vare == true
        nue = 4
        scalee = (nue - 2) / nue * vare
    end

    if my_rank == 0
        if estimate_vara == true
            nub = 4
            scaleb = (nub - 2) / nub * varEffects
        end
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

    xpx_dict = Dict{Int,Vector{Float64}}() # blkID -> vector of xpx for each category
    xArray_dict = Dict{Int, Matrix{Float64}}() # blkID -> TransformedX for each category

    for blk in keys(my_TransformedX_dict)
        Xb = my_TransformedX_dict[blk]
        nMarkerb = size(Xb, 2)
        xpx_dict[blk] = [dot(Xb[:, i], Xb[:, i]) for i in 1:nMarkerb]
        xArray_dict[blk] = Xb
    end

    #output

    betaArray = zeros(my_nsnp)
    alphaArray = zeros(my_nsnp)
    deltaArray = ones(my_nsnp)

    #posterior mean 
    meanAlpha = zeros(my_nsnp)

    if my_rank == 0 
        if estimate_pi == true 
            mean_pi = 0.0
            mean_pi2 = 0.0
        end
        
        if estimate_vara == true 
            # beta effect variance
            mean_vara2 = 0.0
        end

            mean_vara = 0.0

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
        println("Pleiotropy percent: ", pleio_percent)
        println("Heritability for trait 1 (h21): ", h21)
        println("Heritability for trait 2 (h22): ", h22)
        println("Number of ranks: ", nrank)
        println("estimate_vare=$estimate_vare,estimate_vara=$estimate_vara")
        println("estimate_pi=$estimate_pi")
        println("analysis_path=$analysis_path")
        println("data_path=$data_path")
        println("marker effect variance: ", varEffects)
        println("residual variance is: ", round.(vare * nInd, digits=3), "  (R*Ind)")
        println("pi=$startPi")
        println("nMarker=$nMarker, nInd=$nInd, nBlocks=$nBlocks")
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

        # compute total varg for each trait without split into different categories
        totalvarg_blk = fill(0.0, my_nblk) 

        nLoci_vec = 0 # number of causal SNP in this rank
        # sse_vec Computed as dot(\beta, \beta) for whole rank
        sse_vec = 0.0


        #for loop for LD blocks
        for b in 1:my_nblk
            blk = my_blkID[b] #block ID
            xpxc = xpx_dict[blk]
            xArrayc = xArray_dict[blk]
            wArray = my_TransformedY_dict[blk][trait]
            SNPIndexb = my_blkSNPsIndex_dict[blk]
            nEigenb, nMarkerb = size(my_TransformedX_dict[blk]) # q & nsnpb

            vare = vare_blk[b] # VareDn
            invVareDn = 1.0 / vare
            
            invVarEffects = 1 / varEffects
            logVarEffects = log(varEffects)
            
            logPi = log(pi)
            logPiComp = log(1 - pi)
            logDelta0 = logPi

            invLhs = 1.0 / (invVareDn + invVarEffects)
            logInvLhsMsigma = log(invLhs) - logVarEffects

            for marker = 1:nMarkerb
                markerIndex = SNPIndexb[marker] #marker position across my_nsnp
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
                    nLoci_vec += 1
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
            end # end marker loop
            

            my_TransformedY_dict[blk][trait] = wArray

            # compute genetic variance and heritability 
            XAb = xArrayc
            what_array_total = zeros(nEigenb)
            alphaArray_b = alphaArray[SNPIndexb]
            what_array_total = XAb * alphaArray_b
            
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
                pi = samplePi(nLoci_vec_sum, nMarker)
                if iter > burnin
                    mean_pi += (pi - mean_pi) * iIter
                    mean_pi2 += (pi ^ 2 - mean_pi2) * iIter
                end
            end
        end

        # broadcast pi_vec from 0 to other ranks 
        pi = MPI.bcast(pi, 0, comm)
        MPI.Barrier(comm)


        
        ### Step3. sample marker variance ###
        sse_vec = dot(betaArray, betaArray)
        MPI.Barrier(comm)
        # gather the results from all ranks
        sse_vec_sum = MPI.Reduce(sse_vec, +, 0, comm)
        if my_rank == 0
            if estimate_vara == true
                varEffects = (sse_vec_sum + nub * scaleb) / rand(Chisq(nMarker + nub))
            end
            if iter > burnin
                mean_vara += (varEffects - mean_vara) * iIter
                if estimate_vara == true
                    mean_vara2 += (varEffects ^ 2 - mean_vara2) * iIter
                end
            end
        end
        if estimate_vara == true
            #broadcast varEffects_vec from 0 to other ranks
            varEffects = MPI.bcast(varEffects, 0, comm)
            MPI.Barrier(comm)
        end

        # summing genetic variance
        if iter > burnin
            totalvarg = sum(totalvarg_blk) # sum totalvarg_blk to get total genetic variance across all blocks (value)
            totalvarg_sum = MPI.Reduce(totalvarg, +, 0, comm) # sum the results from all ranks

            MPI.Barrier(comm)
            if my_rank == 0
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
                    # println("iter $iter")
                    # println("nQTL: ", nLoci_vec_sum)
                    # println("pi:", pi)
                    
                    writedlm(outfile["pi"], pi, ',')
                    writedlm(outfile["beta_effects_variance"], varEffects, ',')
                    writedlm(outfile["genetic_effects_variance"], totalvarg_sum, ',')
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
