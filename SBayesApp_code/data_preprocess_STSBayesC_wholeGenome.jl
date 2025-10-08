using LinearAlgebra, Distributions, Random, SparseArrays
using DelimitedFiles, DataFrames, CSV, JLD2
using Dates
using InteractiveUtils
using Statistics,ProgressMeter;

data_path = ARGS[1]
nrank     = ARGS[2]
nrank     = parse(Int64, nrank)

map_path = data_path * "LDblock_matched.map"
map = CSV.read(map_path, DataFrame)
@show nrow(map)

# bhat for all SNPs having values
bhat = CSV.read(data_path * "bhatXsqrt2pq_GWASfreq.txt",DataFrame)
bhat[!,:SNP] = String.(bhat[!,:SNP])

blkIDs = unique(map[!,:ldBlockID])
nBlocks = length(blkIDs)

# Create a dictionary of sub-matrices
LDmatrix_dict     = Dict{Int, Matrix{Float64}}()
blkSNPsIndex_dict = Dict{Int, Vector{Int64}}()
bhat_dict         = Dict{Int, Vector{Float64}}()
for blk in blkIDs
    @show blk
    snp_ids_blk = String.(vec(readdlm(data_path * "SNPs_block$blk.csv")))
    nSNPs = length(snp_ids_blk)
    blk_indices = collect(1:nSNPs)
    blkSNPsIndex_dict[blk] = blk_indices 
    LDmatrix_dict[blk]     = readdlm(data_path * "LDmatrix_block$blk.csv");
    # order the values of bhat to let it be consistent with snp_ids_blk 
    indexInBhat = zeros(Int, nSNPs)
    for i in eachindex(snp_ids_blk)
        indexInBhat[i] = findfirst(bhat[!,:SNP] .== snp_ids_blk[i])
    end
    bhat_dict[blk] = bhat[indexInBhat,:bAdj]
end

############################################################
#MPI: split data for different ranks (#rank=#parallel jobs)
############################################################
@show nrank
##############################
# 1. split blkIDs for different ranks
##############################

blkIDs          = sort(blkIDs)
nBlocks         = length(blkIDs)
blkID_each_rank = Vector{Vector{Int}}(undef, nrank)
nPervec         = Int(ceil(nBlocks/nrank))
# Split the vector into sub-vectors
for i in 1:nrank
    start_index = (i-1)*nPervec + 1
    end_index = min(i*nPervec, length(blkIDs))
    blkID_each_rank[i] = blkIDs[start_index:end_index]
end

#save data into data folder
mpi_data_path= data_path * "nrank$nrank/bhatXsqrt2pq_GWASfreq/"
mkpath(mpi_data_path)
for i in 1:nrank
      rankID=i-1 #rank starts from 0, not 1
      writedlm(mpi_data_path*"rank$rankID.blkIDs.txt",blkID_each_rank[i], ',')
end

for i in 1:nrank
    @show i 
    rankID=i-1 #rank starts from 0, not 1
    blkID_ranki=Int.(vec(readdlm(mpi_data_path*"rank$rankID.blkIDs.txt", ',')))  #block ID for this rank
      
    LDmatrix_dict_ranki = Dict(key => LDmatrix_dict[key] for key in blkID_ranki) #get LDmatrix_dict for this rank
    JLD2.save(mpi_data_path*"rank$rankID.LDmatrix_dict.jld2", "my_LDmatrix_dict",LDmatrix_dict_ranki)  #save

    bhat_dict_ranki = Dict(key => bhat_dict[key] for key in blkID_ranki) #get anno_matrix_dict for this rank
    JLD2.save(mpi_data_path*"rank$rankID.bhat_dict.jld2", "my_bhat_dict",bhat_dict_ranki)  #save

    blkSNPsIndex_dict_ranki = Dict(key => blkSNPsIndex_dict[key] for key in blkID_ranki) #get blkSNPsIndex_dict for this rank
    JLD2.save(mpi_data_path*"rank$rankID.blkSNPsIndex_dict.jld2","my_blkSNPsIndex_dict",blkSNPsIndex_dict_ranki)  #save
end


