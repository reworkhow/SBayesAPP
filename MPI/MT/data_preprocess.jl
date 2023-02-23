#helper: get all keys of a dictionary
function get_keys_by_value(dict, value)
      keys_list = []
      for (key, val) in dict
          if val == value
              push!(keys_list, key)
          end
      end
      return keys_list
end


cd("/Users/tianjing/Library/CloudStorage/Box-Box/MTSbayesC/SBayesC_tianjing/parallel")

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
# Preprocess
# sort annotationMat and annotationType, such that continuous annotation is before categorical annotation
# Find the indices of continuous and categorical annotations
is_continuous = (annotationType .== "continue")
is_categorical = (annotationType .== "category")
continuous_indices = findall(is_continuous)
categorical_indices = findall(is_categorical)
    
# Reorder the annotation matrix and type vector
annotationMat = [annotationMat[:, continuous_indices] annotationMat[:, categorical_indices]]
annotationType = [fill("continue", count(is_continuous)) ; fill("category", count(is_categorical))]

# Update the `nCon` and `nCat` variables
nCon = length(continuous_indices)
nCat = length(categorical_indices)

# To convert AnnoIndex to Bool vectors (AnnoIndex for continuous group always = true)
annotationConvertedMat = copy(annotationMat) 
annotationConvertedMat[:, findall(annotationType .== "continue")] .= 1; #continuouse 3.5->true; categorical: 1->true, 0->false
AnnoIndex = [Bool.(annotationConvertedMat[:,c]) for c=1:nCategory]; #this is used to skip sampling marker effects when this marker does not have some annotation

# Split genotypes, annotationMat, annotationConvertedMat to blocks based on id_to_category
anno_df       = DataFrame(annotationMat', names(geno_qc_df));
anno_index_df = DataFrame(annotationConvertedMat', names(geno_qc_df));

geno_data_frames      = Dict()
anno_data_frames      = Dict()
anno_index_data_frame = Dict()
for id in names(geno_qc_df)
    blk = id_to_category[id]   # corresponding blk number 
    if !haskey(geno_data_frames, blk)
        geno_data_frames[blk]     = DataFrame() # genotypes
        anno_data_frames[blk]     = DataFrame() # annotationMat
        anno_index_data_frame[blk]= DataFrame() # annotationConvertedMat
    end
    geno_data_frames[blk][:,id]      = geno_qc_df[:, id]
    anno_data_frames[blk][:,id]      = anno_df[:, id]
    anno_index_data_frame[blk][:,id] = anno_index_df[:, id]
end

# convert dataframe to matrix
geno_matrix_dict       = Dict{Int, Matrix}(name => Matrix(df) for (name, df) in geno_data_frames);
anno_matrix_dict       = Dict{Int, Matrix}(name => Matrix(df)' for (name, df) in anno_data_frames);
anno_index_matrix_dict = Dict{Int, Matrix}(name => Matrix(df)' for (name, df) in anno_index_data_frame);

# Save the SNP index for each blk 
SNPID_all = names(geno_qc_df)
# get SNPs for each blk
blkID = unique(values(id_to_category))
blkSNPsIndex_dict = Dict()
for blk in blkID
    if !haskey(blkSNPsIndex_dict, blk)
        blkSNPsIndex_dict[blk] = findall(in.(SNPID_all, Ref(get_keys_by_value(id_to_category,blk))))
    end
end

############################################################
#MPI: split data for different ranks (#rank=#parallel jobs)
############################################################
nrank=3 #3 parallel jobs

##############################
# 1. split blkID for different ranks
##############################
blkID = sort(blkID)
n_blk_per_rank = length(blkID)÷nrank
split_pos=[1:n_blk_per_rank:length(blkID); length(blkID)+1]
remainder = length(blkID) % nrank # Calculate the remainder of the division
blkID_each_rank = [blkID[(i-1)*n_blk_per_rank+1:i*n_blk_per_rank+(i==nrank)*remainder] for i in 1:nrank]

#save data into data folder
mpi_data_path="/Users/tianjing/Library/CloudStorage/Box-Box/MTSbayesC/SBayesC_tianjing/parallel/nrank$nrank/"
for i in 1:nrank
      rankID=i-1 #rank starts from 0, not 1
      writedlm(mpi_data_path*"rank$rankID.blkID.txt",blkID_each_rank[i], ',')
end

##############################
# 2. split geno_matrix_dict for different ranks
##############################
using JLD2
for i in 1:nrank
      rankID=i-1 #rank starts from 0, not 1
      blkID_ranki=Int.(vec(readdlm(mpi_data_path*"rank$rankID.blkID.txt", ',')))  #block ID for this rank
      
      geno_matrix_dict_ranki = Dict(key => geno_matrix_dict[key] for key in blkID_ranki) #get geno_matrix_dict for this rank
      JLD2.save(mpi_data_path*"rank$rankID.geno_matrix_dict.jld2", "my_geno_matrix_dict",geno_matrix_dict_ranki)  #save

      anno_matrix_dict_ranki = Dict(key => anno_matrix_dict[key] for key in blkID_ranki) #get anno_matrix_dict for this rank
      JLD2.save(mpi_data_path*"rank$rankID.anno_matrix_dict.jld2", "my_anno_matrix_dict",anno_matrix_dict_ranki)  #save

      anno_index_matrix_dict_ranki = Dict(key => anno_index_matrix_dict[key] for key in blkID_ranki) #get anno_index_matrix_dict for this rank
      JLD2.save(mpi_data_path*"rank$rankID.anno_index_matrix_dict.jld2", "my_anno_index_matrix_dict",anno_index_matrix_dict_ranki)  #save

      blkSNPsIndex_dict_ranki = Dict(key => blkSNPsIndex_dict[key] for key in blkID_ranki) #get blkSNPsIndex_dict for this rank
      JLD2.save(mpi_data_path*"rank$rankID.blkSNPsIndex_dict.jld2", "my_blkSNPsIndex_dict",blkSNPsIndex_dict_ranki)  #save
end

#read data:
# my_geno_matrix_dict = JLD2.load(mpi_data_path*"rank$rankID.geno_matrix_dict.jld2")["my_geno_matrix_dict"]



JLD2.save("AnnoIndex.jld2", "AnnoIndex",AnnoIndex)

# AnnoIndex = JLD2.load("AnnoIndex.jld2")["AnnoIndex"]

my_geno_matrix_dict = JLD2.load("nrank3/rank1.geno_matrix_dict.jld2")["my_geno_matrix_dict"]