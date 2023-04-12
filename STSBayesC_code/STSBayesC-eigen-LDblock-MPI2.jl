using MPI
using LinearAlgebra, Distributions, Random, SparseArrays
using DelimitedFiles, DataFrames, CSV, JLD2
using Dates
using InteractiveUtils
using Statistics,ProgressMeter;

code_path = ARGS[1]
folder =  ARGS[2]
data_path = ARGS[3]
nrank = ARGS[4]
nrank = parse(Int64,nrank)
nIter = ARGS[5]
nIter = parse(Int64,nIter)
cd(folder)
include(code_path * "helper.jl")
seed = ARGS[6]
seed = parse(Int64,seed)
@show seed

outFreq=50
estimate_vare=true
estimate_vara=true

vare = 0.5
vary = 1.0
startPi = 0.9

nMarker = 18588
nInd    = 503
nBlocks = 10

function runMPI(;niter=nIter, outFreq=outFreq, seed=seed,
    estimate_vare=estimate_vare, estimate_vara=estimate_vara,
    folder=folder, data_path=data_path,
    vare=vare, vary=vary, pi=startPi, estimatePi=true,
    nMarker=nMarker, n=nInd, nBlocks = nBlocks)

    MPI.Init()
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

    my_TransformedX_dict = JLD2.load(data_path * "nrank$cluster_size.eigen/99Eigen/rank$my_rank.TransformedX_dict.jld2")["my_TransformedX_dict"]
    my_TransformedY_dict = JLD2.load(data_path * "nrank$cluster_size.eigen/99Eigen/rank$my_rank.TransformedY_dict.jld2")["my_TransformedY_dict"]
    my_blkSNPsIndex_dict = JLD2.load(data_path * "nrank$cluster_size.eigen/99Eigen/rank$my_rank.blkSNPsIndex_dict.jld2")["my_blkSNPsIndex_dict"]
    my_blkID = Int.(vec(readdlm(data_path * "nrank$cluster_size.eigen/99Eigen/rank$my_rank.blkIDs.txt", ',')))  #block IDs for this rank
    my_nblk = length(my_blkID)
    my_nsnp = sum(map(length, values(my_blkSNPsIndex_dict)))

    # reorder my_blkSNPsIndex_dict as from 1 to my_nsnp
    # extract smallest values in my_blkSNPsIndex_dict
    my_min = minimum(map(minimum, values(my_blkSNPsIndex_dict)))
    my_blkSNPsIndex_dict = Dict(key => my_blkSNPsIndex_dict[key] .- my_min .+ 1 for key in my_blkID)

    h2 = 0.5
    m = my_nsnp
    varg = vary * h2
    vare = vary * (1 - h2)
    varEffects = varg / (nMarker * (1 - pi))

    α = zeros(m) # a vector of estimated marker effects
    β = zeros(m)
    δ = ones(m)

    α_mcmc = zeros(niter, m)
    nub = 4
    nue = 4
    scaleb = (nub - 2) / nub * varEffects
    scalee = (nue - 2) / nue * vare

    if my_rank == 0 # only root rank save the output
        nOutput = Int(floor(niter / outFreq))
        pi_mcmc = zeros(niter)
        varEffects_mcmc = zeros(niter)
        vare_mcmc = zeros(niter)
        hsq_mcmc = zeros(niter)
        keptIter = zeros(nOutput, 6)
    end

    my_Wcorr_dict = deepcopy(my_TransformedY_dict)
    vare_blk = zeros(my_nblk)
    varg_blk = zeros(my_nblk)
    iout = 1

    ############################################################################
	# print information
	############################################################################
	if my_rank == 0
		println("---------------- Summary Start --------------")
		println("nIter=$niter, outFreq=$outFreq, seed=$seed")
        println("estimate_vare=$estimate_vare,estimate_vara=$estimate_vara,estimatePi=$estimatePi")
        println("folder=$folder,data_path=$data_path")
        println("vare=$vare,varEffects=$varEffects, pi=$pi")
        println("nMarker=$nMarker, nInd=$n, nBlocks=$nBlocks")
		time_start = now()
		println("Start time: ", time_start)
		println("---------------- Summary End ----------------")
	end
	MPI.Barrier(comm)
    println("In rank$my_rank, there are $my_nblk LD blocks, and $my_nsnp SNPs in total.")
	MPI.Barrier(comm)

    @showprogress "running MCMC ..." for iter in 1:niter
        nsnp_till_last_blk = 0
        # sampling SNP effects
        logPi = log(pi)
        logPiComp = log(1 - pi)
        invVarEffects = 1 / varEffects
        logVarEffects = log(varEffects)
        logDelta0 = logPi

        nloci = 0
        αTα   = 0.0

        for b in 1:my_nblk
            blk = my_blkID[b]
            Q = my_TransformedX_dict[blk]
            wcorr = my_Wcorr_dict[blk]
            blkSNPsIndex = my_blkSNPsIndex_dict[blk]
            nsnpb = length(blkSNPsIndex)
            q = length(wcorr)  # number of eigen vectors

            # initial what
            what = zeros(q)

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
                    BLAS.axpy!(oldAlpha - α[i+nsnp_till_last_blk], x, wcorr)
                    # update what 
                    what += α[i+nsnp_till_last_blk] * x
                    #update nloci 
                    nloci += 1
                    #update αTα
                    αTα += α[i+nsnp_till_last_blk]^2
                else
                    if (oldAlpha != 0.0)
                        BLAS.axpy!(oldAlpha, x, wcorr)
                    end
                    δ[i+nsnp_till_last_blk] = 0
                    β[i+nsnp_till_last_blk] = randn() * sqrt(varEffects)
                    α[i+nsnp_till_last_blk] = 0
                end
            end   # end of marker loop 
            nsnp_till_last_blk += nsnpb
            my_Wcorr_dict[blk] = wcorr
            vare_blk[b] = (n * dot(wcorr, wcorr) + nue * scalee) / rand(Chisq(q + nue))

            # compute genetic variance and heritability
            varg_blk[b] = dot(what, what)      # y = Xβ + e 
        end # end of blk loop

        α_mcmc[iter, :] = α

        # summing residual variance
        vare_blk_sum = sum(vare_blk)

        # summing genetic variance
        varg = sum(varg_blk)     

        # gather the results
        MPI.Barrier(comm)
        nonzero_gather = MPI.Gather(nloci, 1, 0, comm)
        αTα_gather = MPI.Gather(αTα, 1, 0, comm)
        vare_blk_sum_gather = MPI.Gather(vare_blk_sum, 1, 0, comm)
        varg_gather = MPI.Gather(varg, 1, 0, comm)
        MPI.Barrier(comm)

        # sampling pi
        if (estimatePi)
            if my_rank == 0
                nloci_sum = sum(nonzero_gather)
                pi = samplePi(nloci_sum, nMarker)
                pi_mcmc[iter] = pi
            end
        end

        # sampling SNP effect variance
        if my_rank == 0
            αTα_sum = sum(αTα_gather)
            varEffects = (αTα_sum + nub * scaleb) / rand(Chisq(nloci_sum + nub))
            varEffects_mcmc[iter] = varEffects
        end

        # sampling residual variance
        if my_rank ==0 
            vare_blk_sum_sum = sum(vare_blk_sum_gather)
            vare = vare_blk_sum_sum/nBlocks
            vare_mcmc[iter] = vare
        end

        # compute genetic variance and heritability 
        if my_rank == 0
            varg_sum = sum(varg_gather)
            hsq = varg_sum / vary
            hsq_mcmc[iter] = hsq
        end
        MPI.Barrier(comm)

        #broadcast from rank 0 to other ranks
        varEffects = MPI.bcast(varEffects, 0, comm) 
        pi         = MPI.bcast(pi, 0, comm)
        if estimate_vare == true
            vare   = MPI.bcast(vare, 0, comm)
        end
        MPI.Barrier(comm)

        if my_rank == 0 
            if  iter % outFreq == 0
                keptIter[iout, :] = [pi, nloci_sum, varEffects, hsq, vare, varg_sum]
                println("iter $iter, pi = $pi, nloci = $nloci_sum, varEffects = $varEffects, hsq = $hsq, vare = $vare, varg = $varg_sum.")
                iout += 1
            end
        end
        MPI.Barrier(comm)
    end # end of mcmc iter loop

    if my_rank == 0 
        df_name = ["pi", "nloci", "varEffects", "hsq", "Vare", "Varg"]
        par = mean(keptIter, dims=1)
        posterior_means = DataFrame(:parameter => df_name, :mean => vec(par))
        CSV.write(folder * "posterior_means_bcpi.csv", posterior_means)

        # save mcmc samples
        writedlm(folder * "mcmcPi.txt",pi_mcmc)
        writedlm(folder * "mcmcVarEffects.txt",varEffects_mcmc)
        writedlm(folder * "mcmcVare.txt",vare_mcmc)
        writedlm(folder * "mcmcHsq.txt",hsq_mcmc)
    end

    writedlm(folder * "mcmc_Alpha.rank$my_rank.txt",α_mcmc)
    writedlm(folder * "mean_Alpha.rank$my_rank.txt",vec(mean(α_mcmc, dims=1)))

    ############################################################################
	# Step5. show running time
	############################################################################
	if my_rank ==0
		time_end = now()
		time_diff = (time_end-time_start).value/60000 #milliseconds to min
		println("End time: ", time_end)
		println("Running Time (min): ", time_diff)
    end
    MPI.Finalize()
end #end of function


#run MPI
runMPI()


        
    