cd("/Users/tianjing/Library/CloudStorage/Box-Box/MTSbayesC/SBayesC_tianjing/parallel")

include("helper.jl")

using Random,Statistics,LinearAlgebra,Plots,Distributions,DataFrames,ProgressMeter,DelimitedFiles,CSV;

# Input data 
#read genotypes
data_path = ""
geno_qc_df = CSV.read(data_path*"geno_n503_p18588_realSnpID.QC.scale.csv",DataFrame) # genotypes (after QC)
geno_qc_df = geno_qc_df[:,2:end]; #1st column is ID
nInd,nMarker = size(geno_qc_df);
geno_qc_df

#simulate annotation data
a1 = sample([1.,0.],nMarker)
a2 = sample([2.,2.25,2.5,2.75,3.],nMarker)
annotationMat = [a1 a2] # matrix input
annotationType = ["category","continue"];
nCategory = length(annotationType);

#read map file (last column is block ID)
map_chr22=CSV.read(data_path*"g1000_eur.chr22.hapmap3.newblk.map",DataFrame)
# split into each blocks
id_to_category = Dict(zip(map_chr22[:, :snpID], map_chr22[:, :newBlk]));


#################################### HELPER ####################################
# # Preprocess
# # sort annotationMat and annotationType, such that continuous annotation is before categorical annotation
# # Find the indices of continuous and categorical annotations
# is_continuous = (annotationType .== "continue")
# is_categorical = (annotationType .== "category")
# continuous_indices = findall(is_continuous)
# categorical_indices = findall(is_categorical)
    
# # Reorder the annotation matrix and type vector
# annotationMat = [annotationMat[:, continuous_indices] annotationMat[:, categorical_indices]]
# annotationType = [fill("continue", count(is_continuous)) ; fill("category", count(is_categorical))]

# # Update the `nCon` and `nCat` variables
# nCon = length(continuous_indices)
# nCat = length(categorical_indices)

# # To convert AnnoIndex to Bool vectors (AnnoIndex for continuous group always = true)
# annotationConvertedMat = copy(annotationMat) 
# annotationConvertedMat[:, findall(annotationType .== "continue")] .= 1; #continuouse 3.5->true; categorical: 1->true, 0->false
# AnnoIndex = [Bool.(annotationConvertedMat[:,c]) for c=1:nCategory]; #this is used to skip sampling marker effects when this marker does not have some annotation

# # Split genotypes, annotationMat, annotationConvertedMat to blocks based on id_to_category
# anno_df       = DataFrame(annotationMat', names(geno_qc_df));
# anno_index_df = DataFrame(annotationConvertedMat', names(geno_qc_df));

# geno_data_frames      = Dict()
# anno_data_frames      = Dict()
# anno_index_data_frame = Dict()
# for id in names(geno_qc_df)
#     blk = id_to_category[id]   # corresponding blk number 
#     if !haskey(geno_data_frames, blk)
#         geno_data_frames[blk]     = DataFrame() # genotypes
#         anno_data_frames[blk]     = DataFrame() # annotationMat
#         anno_index_data_frame[blk]= DataFrame() # annotationConvertedMat
#     end
#     geno_data_frames[blk][:,id]      = geno_qc_df[:, id]
#     anno_data_frames[blk][:,id]      = anno_df[:, id]
#     anno_index_data_frame[blk][:,id] = anno_index_df[:, id]
# end

# # convert dataframe to matrix
# geno_matrix_dict       = Dict{Int, Matrix}(name => Matrix(df) for (name, df) in geno_data_frames);
# anno_matrix_dict       = Dict{Int, Matrix}(name => Matrix(df)' for (name, df) in anno_data_frames);
# anno_index_matrix_dict = Dict{Int, Matrix}(name => Matrix(df)' for (name, df) in anno_index_data_frame);

# # Save the SNP index for each blk 
# SNPID_all = names(geno_qc_df)
# # get SNPs for each blk
# blkID = unique(values(id_to_category))
# blkSNPsIndex_dict = Dict()
# for blk in blkID
#     if !haskey(blkSNPsIndex_dict, blk)
#         blkSNPsIndex_dict[blk] = findall(in.(SNPID_all, Ref(get_keys_by_value(id_to_category,blk))))
#     end
# end
# ############################################################################################################


# Input data:
    # y_all (all), geno_matrix_dict (rank), blkSNPsIndex_dict (rank), anno_matrix_dict (rank), anno_index_matrix_dict (rank), annotationType (all),
    # blkID (rank), nCon (all), nCat (all), AnnoIndex,
# MCMC parameters:
    # Rprior, Gprior, startPi, 
    # nIter, outFreq, seed, estimate_vare, estimate_varg, folder
    # nMarker,nCategory,nInd



y1     = vec(readdlm(data_path * "y1.txt"))
y2     = vec(readdlm(data_path * "y2.txt"))
y_all = [y1, y2];

Rprior = [0.3 0.23; 0.23 0.72]
Gprior = [0.0038 0.0015; 0.0015 0.0016]     # marker effect covariance
startPi=Dict([1.0; 1.0]=>0.7,[1.0; 0.0]=>0.1,[0.0; 1.0]=>0.1,[0.0; 0.0]=>0.1)

nIter=10
outFreq=1
seed=123
estimate_vare=false
estimate_varg=true
folder="";

Random.seed!(seed)

#Initialization
G_vec = [Gprior for c in 1:nCategory]
Rinv  = inv(Rprior)
Pi    = [startPi for c in 1:nCategory]

nTraits= length(y_all)

# hyper-parameters for G and R
df_R         = 4 + nTraits
df_G         = 4 + nTraits
scale_R      = Rprior * (df_R - nTraits - 1)
scale_G      = Gprior/nMarker * (df_G - nTraits - 1)

wArray    = deepcopy(y_all)
iout      = 1  

betaArray = [zeros(nMarker*nCategory) for t in 1:nTraits] #-> ordered by annotation groups (##better to be vector of vector of vector?)
alphaArray= [zeros(nMarker*nCategory) for t in 1:nTraits]
deltaArray= [ones(nMarker*nCategory) for t in 1:nTraits]
meanAlpha = [zeros(nMarker*nCategory) for t in 1:nTraits]

if estimate_varg == true
	meanG     = [zeros(nTraits,nTraits) for c in 1:nCategory]
end

nOutput        = Int(floor(nIter/outFreq))
mcmc_alphaArray= [zeros(nMarker*nCategory,nOutput) for t in 1:nTraits]
mcmc_meanAlpha = [zeros(nMarker*nCategory,nOutput) for t in 1:nTraits]
@showprogress "running MCMC ..." for iter = 1:nIter
    iIter = 1.0/iter
    #for loop for LD blocks
    for blk in blkID
        ########################### move out of iteration loop? ###########################
        Xb                      = geno_matrix_dict[blk]  #genotype for current blk, n-by-p (p is #snp in this blk)
        annotationMatb          = anno_matrix_dict[blk]  #annotation data for current blk, p-by-C
        annotationConvertedMatb = anno_index_matrix_dict[blk]  #annotation true/false for current blk, p-by-C
        SNPIndexb               = blkSNPsIndex_dict[blk] #snp pos for current blk, p-by-1

        AnnoIndexb = [Bool.(annotationConvertedMatb[:,c]) for c=1:nCategory] #annotation true/false, vector of length C, each element is vector of length p
        aArrayb    = [annotationMatb[:,c] for c=1:nCategory] #annotation data, vector of length C, vector of length C, each element is vector of length p
        nMarkerb  = length(SNPIndexb)

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
        xArray_vec[end] = [Xb[:,i] for i in 1:nMarkerb]
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
            annoindexc = AnnoIndexb[cat]
            if cat <= nCon+1 # continuous group + categorical group (after nCon+1, xArrayc/xpxc will not change)
                xArrayc = xArray_vec[cat]
                xpxc    = xpx_vec[cat]
            end
            for marker=1:nMarkerb
                true_marker_num = SNPIndexb[marker] #marker position
                if annoindexc[marker] # skip sampling if the marker is not in the category
                    markerIndex = (cat-1)*nMarker + true_marker_num # exact position in betaArray; betaArray[trait1]: nMarker*nCategory-by-1
                    x = xArrayc[marker]
                    for trait = 1:nTraits
                            β[trait]= betaArray[trait][markerIndex]
                        oldα[trait] = newα[trait] = alphaArray[trait][markerIndex]
                           δ[trait] = deltaArray[trait][markerIndex]
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
                        if(rand()<probDelta1) #δj=1
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
                        betaArray[trait][markerIndex]  = β[trait] 
                        deltaArray[trait][markerIndex] = δ[trait]
                        alphaArray[trait][markerIndex] = newα[trait]
                        meanAlpha[trait][markerIndex] += (newα[trait]-meanAlpha[trait][markerIndex])*iIter
                    end
                end # end if loop  
            end # end marker loop
        end  # end annotation loop
    end # end blk loop 
    
    ### Step2. sample Pi ###
    for cat=1:nCategory
        deltaArrayc = [deltaArray[t][((cat-1)*nMarker + 1):((cat-1)*nMarker + nMarker)][AnnoIndex[cat]] for t in 1:nTraits]
        Pi[cat] = mysamplePi(deltaArrayc,Pi[cat]) #annotation specific pi
    end
    
    ### Step3. sample variance ###
    if estimate_varg == true
        for cat=1:nCategory
            betaArrayc = [betaArray[t][((cat-1)*nMarker + 1):((cat-1)*nMarker + nMarker)][AnnoIndex[cat]] for t in 1:nTraits]
            G_vec[cat]=sample_variance(betaArrayc, nMarker, df_G, scale_G, false, false) #annotation specific G
        end
    end
    
    if iter%outFreq ==0
        # marker effects
        for trait = 1:nTraits
            mcmc_alphaArray[trait][:,iout] = alphaArray[trait]
            mcmc_meanAlpha[trait][:,iout]  = meanAlpha[trait]
        end
        iout +=1
	end
    
    if estimate_varg == true
        for cat in 1:nCategory
            meanG[cat] += (G_vec[cat] - meanG[cat])*iIter
        end
	end
end # end iteration loop

if estimate_varg == true
    for cat in 1:nCategory
        writedlm(folder * "estG"*string(cat)*".txt",meanG[cat])
    end
end

writedlm(folder * "mcmc_meanAlpha1.txt",mcmc_meanAlpha[1])
writedlm(folder * "mcmc_meanAlpha2.txt",mcmc_meanAlpha[2])
writedlm(folder * "mcmc_alpha1.txt",mcmc_alphaArray[1])
writedlm(folder * "mcmc_alpha2.txt",mcmc_alphaArray[2])