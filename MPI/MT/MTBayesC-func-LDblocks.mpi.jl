
# Input data:
    # y_all (all), geno_matrix_dict (rank), blkSNPsIndex_dict (rank), anno_matrix_dict (rank), anno_index_matrix_dict (rank), annotationType (all),
    # blkID (rank), nCon (all), nCat (all), AnnoIndex,
# MCMC parameters:
    # Rprior, Gprior, startPi, 
    # nIter, outFreq, seed, estimate_vare, estimate_varg, folder
    # nMarker,nCategory,nInd


using MPI
using LinearAlgebra, Distributions, Random, SparseArrays
using DelimitedFiles, DataFrames, CSV, JLD2
using Dates
using InteractiveUtils
using Random,Statistics,LinearAlgebra,Plots,Distributions,DataFrames,ProgressMeter,DelimitedFiles,CSV;

include("helper.jl")

cd("/group/qtlchenggrp/tianjing/MTBayesC/parallel")

nIter=10
outFreq=1
seed=123
estimate_vare=false
estimate_varg=true
folder=""
data_path="/group/qtlchenggrp/tianjing/MTBayesC/parallel/"
nCon=1
nCategory=1
nMarker=18588

Rprior = [0.3 0.23; 0.23 0.72]
Gprior = [0.0038 0.0015; 0.0015 0.0016]     # marker effect covariance
startPi=Dict([1.0; 1.0]=>0.7,[1.0; 0.0]=>0.1,[0.0; 1.0]=>0.1,[0.0; 0.0]=>0.1)

nCategory=2

function t(;nIter=nIter,outFreq=outFreq,seed=seed,
            estimate_vare=estimate_vare,estimate_varg=estimate_varg,
            folder=folder,data_path=data_path,
            Rprior=Rprior,Gprior=Gprior,startPi=startPi,
            nCategory=nCategory,nCon=nCon,nMarker=nMarker)

      MPI.Init()
      comm         = MPI.COMM_WORLD
      my_rank      = MPI.Comm_rank(comm) #current rank, e.g., 0/1/2/3, root=0
      cluster_size = MPI.Comm_size(comm) #number of all processes/rank, e.g., 4
      @show my_rank,cluster_size
      MPI.Barrier(comm)

      # set seed in different rank
	Random.seed!(seed+my_rank)
	MPI.Barrier(comm)

      ############################################################################
      #read data in current rank
      ############################################################################
      my_geno_matrix_dict = JLD2.load(data_path * "nrank$cluster_size/rank$my_rank.geno_matrix_dict.jld2")["my_geno_matrix_dict"]
      my_anno_matrix_dict = JLD2.load(data_path * "nrank$cluster_size/rank$my_rank.anno_matrix_dict.jld2")["my_anno_matrix_dict"]
      my_anno_index_matrix_dict = JLD2.load(data_path * "nrank$cluster_size/rank$my_rank.anno_index_matrix_dict.jld2")["my_anno_index_matrix_dict"]
      my_blkSNPsIndex_dict = JLD2.load(data_path * "nrank$cluster_size/rank$my_rank.blkSNPsIndex_dict.jld2")["my_blkSNPsIndex_dict"]
      my_blkID = Int.(vec(readdlm(data_path*"nrank$cluster_size/rank$my_rank.blkID.txt", ',')))  #block ID for this rank
      my_nblk = length(my_blkID)
      my_nsnp = sum(map(length, values(my_blkSNPsIndex_dict))) #total number of SNPs in this rank (=sum of #snp in all blocks of this rank)

      y1      = vec(readdlm(data_path * "y1.txt"))
      y2      = vec(readdlm(data_path * "y2.txt"))
      y_all   = [y1, y2]; #vector of vector

      nsnp_all_rank = MPI.Gather(my_nsnp, 0, comm) #collect number of SNPs in each rank, e.g., [1000,2000,1500,...]

      if my_rank==0 #AnnoIndex of all markers
            AnnoIndex = JLD2.load("AnnoIndex.jld2")["AnnoIndex"] #vector of length #category, each element is a vector of length #allSNPs
      end

      ############################################################################
      # initialization in all ranks
      ############################################################################
      nTraits = length(y_all)
      wArray  = deepcopy(y_all) #assume each rank are independent, so they will have different ycorr/wArray. Should we have different wArray for each block?

      G_vec = [Gprior for _ in 1:nCategory]
      Rinv  = inv(Rprior)  #fixed
      Pi    = [startPi for _ in 1:nCategory]
      
      my_betaArray = [zeros(my_nsnp,nCategory) for _ in 1:nTraits] #vector of matrix, each matrix is #SnpInThisRank-by-#annotation
      my_alphaArray= [zeros(my_nsnp,nCategory) for _ in 1:nTraits] # e.g., for rank0: #SnpInThisRank = #snp_blk1 + #snp_blk2...
      my_deltaArray= [ones(my_nsnp,nCategory)  for _ in 1:nTraits]
      my_meanAlpha = [zeros(my_nsnp,nCategory) for _ in 1:nTraits]

      betaArray_vbuf  = [] #vector of VBuffer
      alphaArray_vbuf = [] #vector of VBuffer
      deltaArray_vbuf = [] #vector of VBuffer
      meanAlpha_vbuf  = [] #vector of VBuffer

      ############################################################################
	# initial data and buffer in rank 0
	############################################################################
	if my_rank == 0
		counts = nsnp_all_rank*nCategory  # number of elements in each rank,e.g., [1000*3,2000*3,1500*3,...]

            betaArray  = [zeros(nMarker,nCategory) for _ in 1:nTraits]
            alphaArray = [zeros(nMarker,nCategory) for _ in 1:nTraits]
            deltaArray = [ones(nMarker,nCategory) for _ in 1:nTraits]
            meanAlpha = [zeros(nMarker,nCategory) for _ in 1:nTraits]

            for tt in 1:nTraits
                  push!(betaArray_vbuf,  VBuffer(betaArray[tt],  counts) )
                  push!(alphaArray_vbuf, VBuffer(alphaArray[tt], counts) )
                  push!(deltaArray_vbuf, VBuffer(deltaArray[tt], counts) )
                  push!(meanAlpha_vbuf,  VBuffer(meanAlpha[tt],  counts) )
            end

            #G
            df_G    = 4 + nTraits
            scale_G = Gprior/nMarker * (df_G - nTraits - 1)
            if estimate_varg == true
                  meanG     = [zeros(nTraits,nTraits) for _ in 1:nCategory]
            end

            outfile1 = open(folder * "mcmc_meanAlpha_trait1.txt", "w")
		outfile2 = open(folder * "mcmc_meanAlpha_trait2.txt", "w")
            outfile3 = open(folder * "mcmc_alpha_trait1.txt", "w")
            outfile4 = open(folder * "mcmc_alpha_trait2.txt", "w")
      else
            for _ in 1:nTraits
                  push!(betaArray_vbuf,  VBuffer(nothing) )
                  push!(alphaArray_vbuf, VBuffer(nothing) )
                  push!(deltaArray_vbuf, VBuffer(nothing) )
                  push!(meanAlpha_vbuf,  VBuffer(nothing) )
            end
      end

      ############################################################################
	# print information
	############################################################################
	if my_rank == 0
		println("---------------- Summary Start --------------")
		println("nIter=$nIter, outFreq=$outFreq, seed=$seed")
            println("estimate_vare=$estimate_vare,estimate_varg=$estimate_varg")
            println("folder=$folder,data_path=$data_path")
            println("Rprior=$Rprior,Gprior=$Gprior, startPi=$startPi")
            println("nCategory=$nCategory,nCon=$nCon,nMarker=$nMarker")
		t_start = now()
		println("Start time: ", t_start)
		println("---------------- Summary End ----------------")
	end
	MPI.Barrier(comm)
      println("In rank$my_rank, there are $my_nblk LD blocks, and $my_nsnp SNPs in total.")
	MPI.Barrier(comm)


      ############################################################################
	# MCMC iterations start (run in current rank)
	############################################################################
      for iter = 1:nIter
            iIter = 1.0/iter
            nsnp_till_last_blk=0 #within this rank
            #for loop for LD blocks
            for blk in my_blkID
                  ########################### move out of iteration loop? ###########################
                  Xb                      = my_geno_matrix_dict[blk]       #genotype for current blk, n-by-p (p is #snp in this blk)
                  annotationMatb          = my_anno_matrix_dict[blk]       #annotation data for current blk, p-by-C
                  annotationConvertedMatb = my_anno_index_matrix_dict[blk] #annotation true/false for current blk, p-by-C
                  SNPIndexb               = my_blkSNPsIndex_dict[blk]      #true snp index for current blk, p-by-1
                  
                  AnnoIndexb = [Bool.(annotationConvertedMatb[:,c]) for c=1:nCategory] #annotation true/false, vector of length C, each element is vector of length p
                  aArrayb    = [annotationMatb[:,c] for c=1:nCategory] #annotation data, vector of length C, each element is vector of length p
                  nMarkerb   = length(SNPIndexb) # #snp in this block

                  # save xpx for continue and category annotation
                  # the last element is for category annotation
                  xpx_vec    = Vector{Vector{Float64}}(undef, nCon+1)
                  xArray_vec = Vector{Vector{Array{Float64,1}}}(undef, nCon+1)
            
                  for c in 1:nCon # xpx and xArray will change for continuous group 
                        xpx_vec[c]    = [aArrayb[c][i]^2 * dot(Xb[:,i], Xb[:,i]) for i in 1:nMarkerb]
                        xArray_vec[c] = [aArrayb[c][i]*Xb[:,i] for i in 1:nMarkerb]
                  end

                  # xpx/xArray will not change for categorical group if the marker effects of SNPs not in the group = 0
                  xpx_vec[end]    = [dot(Xb[:,i], Xb[:,i]) for i in 1:nMarkerb] 
                  xArray_vec[end] = [Xb[:,i] for i in 1:nMarkerb]   #future improve: use reference, instead of copy 
                  ########################### move out of iteration loop? ########################### 
                  
                  β         = zeros(nTraits)  
                  newα      = zeros(nTraits)  
                  oldα      = zeros(nTraits)  
                  w         = zeros(nTraits)  
                  δ         = ones(nTraits)  
            
                  # initialize xArrayc/xpxc 
                  xArrayc = similar(xArray_vec[1])
                  xpxc    = similar(xpx_vec[1])
                  
                  for cat=1:nCategory
                        annoindexc = AnnoIndexb[cat]  #annotation true/false for this category, vector of length p
                        if cat <= nCon+1 # continuous group + categorical group (after nCon+1, xArrayc/xpxc will not change)
                              xArrayc = xArray_vec[cat] #x for this category, vector of length p
                              xpxc    = xpx_vec[cat]    #x'x for this category, vector of length p
                        end
                        for marker=1:nMarkerb
                              # true_marker_num = SNPIndexb[marker] #true marker position
                              if annoindexc[marker] # skip sampling if the marker is not in the category
                                    # markerIndex = (cat-1)*nMarker + true_marker_num # exact position in betaArray; betaArray[trait1]: nMarker*nCategory-by-1
                                    x = xArrayc[marker]
                                    for trait = 1:nTraits
                                          β[trait]= my_betaArray[trait][nsnp_till_last_blk+marker,cat]
                                          oldα[trait] = newα[trait] = my_alphaArray[trait][nsnp_till_last_blk+marker,cat]
                                          δ[trait] = my_deltaArray[trait][nsnp_till_last_blk+marker,cat]
                                          w[trait] = dot(x,wArray[trait])+xpxc[marker]*oldα[trait] #w=xj'(ycorr+xj*αj) scaler
                                    end
                                    Ginv = inv(G_vec[cat])
                                    for k=1:nTraits
                                          Ginv11 = Ginv[k,k]
                                          nok    = deleteat!(collect(1:nTraits),k)
                                          Ginv12 = Ginv[k,nok] #this is not row vector!!, so this is Ginv21
                                          C11    = Ginv11+Rinv[k,k]*xpxc[marker]
                                          C12    = Ginv12+xpxc[marker]*Diagonal(δ[nok])*Rinv[k,nok] #C21
                                          #when δj=0
                                          invLhs0  = 1/Ginv11
                                          rhs0     = - dot(Ginv12,β[nok])
                                          gHat0    = rhs0*invLhs0
                                          #when δj=1
                                          invLhs1  = 1/C11
                                          rhs1     = w'*Rinv[:,k]-C12'β[nok]  #here the w' is in paper: w'm=xj'(ycorr+xj*βj)
                                          gHat1    = rhs1*invLhs1

                                          d0    = copy(δ)
                                          d1    = copy(δ)
                                          d0[k] = 0.0
                                          d1[k] = 1.0

                                          #sample δj
                                          BigPi = Pi[cat]
                                          logDelta0  = -0.5*(log(Ginv11)- gHat0^2*Ginv11) + log(BigPi[d0]) #logPi
                                          logDelta1  = -0.5*(log(C11)-gHat1^2*C11) + log(BigPi[d1]) #logPiComp
                                          probDelta1 =  1.0/(1.0+exp(logDelta0-logDelta1))
                                          
                                          #sample marker effects
                                          if (rand()<probDelta1) #δj=1
                                                δ[k] = 1
                                                β[k] = newα[k] = gHat1 + randn()*sqrt(invLhs1)
                                                wArray[k] = wArray[k] + x*(oldα[k]-newα[k])
                                          else
                                                β[k] = gHat0 + randn()*sqrt(invLhs0)
                                                δ[k] = 0
                                                newα[k] = 0
                                                if oldα[k] != 0
                                                      wArray[k] = wArray[k] + x*oldα[k] #newα[k]=0
                                                end
                                          end
                                    end
                                    
                                    for trait = 1:nTraits
                                          my_betaArray[trait][nsnp_till_last_blk+marker,cat]  = β[trait] 
                                          my_deltaArray[trait][nsnp_till_last_blk+marker,cat] = δ[trait]
                                          my_alphaArray[trait][nsnp_till_last_blk+marker,cat] = newα[trait]
                                          my_meanAlpha[trait][nsnp_till_last_blk+marker,cat] += (newα[trait]-my_meanAlpha[trait][nsnp_till_last_blk+marker,cat])*iIter
                                    end
                              end # end if loop  
                        end # end marker loop
                  end  # end annotation loop

                  nsnp_till_last_blk += nMarkerb
            end # end blk loop

            #send my_betaArray, my_deltaArray, my_alphaArray, my_meanAlpha to rank0
            for trait =1:nTraits
                  MPI.Gatherv!(my_betaArray[trait],  betaArray_vbuf[trait],  0, comm)
                  MPI.Gatherv!(my_deltaArray[trait], deltaArray_vbuf[trait], 0, comm)
                  MPI.Gatherv!(my_alphaArray[trait], alphaArray_vbuf[trait], 0, comm)
                  MPI.Gatherv!(my_meanAlpha[trait],  meanAlpha_vbuf[trait],  0, comm)
            end
		MPI.Barrier(comm)
            
            ########################################################################
	      # Step3. sample Pi, G_vec in rank 0
		########################################################################
            if my_rank==0
                  ### change order of betaArray, deltaArray, alphaArray, meanAlpha
                  #this is because MPI fill the element by row, but we hope to stacking matrix vertically
                  for t in 1:nTraits
                        betaArray[t]  = Matrix(reshape(vec(betaArray[t]),nCategory,:)')   #change to stacking matrix vertically
                        deltaArray[t] = Matrix(reshape(vec(deltaArray[t]),nCategory,:)')  #change to stacking matrix vertically
                        alphaArray[t] = Matrix(reshape(vec(alphaArray[t]),nCategory,:)')  #change to stacking matrix vertically
                        meanAlpha[t] = Matrix(reshape(vec(meanAlpha[t]),nCategory,:)')    #change to stacking matrix vertically
                  end

                  ### sample Pi
                  for cat=1:nCategory
                        deltaArrayc = [deltaArray[t][:,cat][AnnoIndex[cat]] for t in 1:nTraits]
                        Pi[cat] = mysamplePi(deltaArrayc,Pi[cat]) #annotation specific pi
                  end
            
                  ### sample variance 
                  if estimate_varg == true
                        for cat=1:nCategory
                              betaArrayc = [betaArray[t][:,cat][AnnoIndex[cat]] for t in 1:nTraits]
                              G_vec[cat]=sample_variance(betaArrayc, nMarker, df_G, scale_G, false, false) #annotation specific G
                              meanG[cat] += (G_vec[cat] - meanG[cat])*iIter #update
                        end
                  end
            end
            MPI.Barrier(comm)

            #broadcast from rank0 to other ranks
            G_vec = MPI.bcast(G_vec, 0, comm) 
            Pi    = MPI.bcast(Pi, 0, comm) 

            ########################################################################
		# Step4. save mcmc samples in rank 0
		########################################################################
            if iter%outFreq==0 
                  if my_rank==0
                        # marker effects      
                        writedlm(outfile1, alphaArray[1]'  ,  ',')
				writedlm(outfile2, alphaArray[2]'  ,  ',') 
				writedlm(outfile3, meanAlpha[1]' ,  ',') 
				writedlm(outfile4, meanAlpha[2]' ,  ',') 
                  end
            end
            MPI.Barrier(comm) #new added

      end # end MCMC iteration loop


      ############################################################################
	# Step5. show running time
	############################################################################
	if my_rank ==0
		t_end = now()
		t_diff = (t_end-t_start).value/60000 #milliseconds to min
		println("End time: ", t_end)
		println("Running Time (min): ", t_diff)
	end

      #save other results
      if my_rank==0
            if estimate_varg == true
                  for cat in 1:nCategory
                        writedlm(folder * "estG"*string(cat)*".txt",meanG[cat])
                  end
            end
      end


      MPI.Finalize()
end #end t()


#run MPI
t()