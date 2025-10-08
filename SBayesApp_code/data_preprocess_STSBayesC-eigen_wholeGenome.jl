using LinearAlgebra, Distributions, Random, SparseArrays
using DelimitedFiles, DataFrames, CSV, JLD2
using Dates
using InteractiveUtils
using Statistics,ProgressMeter;


data_path = ARGS[1]
nrank     = ARGS[2]
nrank     = parse(Int64, nrank)
@show nrank

mpi_data_path= data_path * "nrank$nrank.eigen/bhatXsqrt2pq_GWASfreq/995Eigen/"
mkpath(mpi_data_path)
for i in 1:nrank
    @show i 
    rankID=i-1 #rank starts from 0, not 1
    blkID_ranki             =Int.(vec(readdlm(data_path * "nrank$nrank/bhatXsqrt2pq_GWASfreq/rank$rankID.blkIDs.txt", ',')))  
    writedlm(mpi_data_path*"rank$rankID.blkIDs.txt",blkID_ranki, ',')
    
    LDmatrix_dict_ranki     = JLD2.load(data_path * "nrank$nrank/bhatXsqrt2pq_GWASfreq/rank$rankID.LDmatrix_dict.jld2")["my_LDmatrix_dict"]
    bhat_dict_ranki         = JLD2.load(data_path * "nrank$nrank/bhatXsqrt2pq_GWASfreq/rank$rankID.bhat_dict.jld2")["my_bhat_dict"]
    blkSNPsIndex_dict_ranki = JLD2.load(data_path * "nrank$nrank/bhatXsqrt2pq_GWASfreq/rank$rankID.blkSNPsIndex_dict.jld2")["my_blkSNPsIndex_dict"]
    TransformedX_dict_ranki = Dict{Int, Matrix{Float64}}()
    TransformedY_dict_ranki = Dict{Int, Vector{Float64}}()
    for blk in blkID_ranki
        @show blk
        LDmatrix_blk               = LDmatrix_dict_ranki[blk]
        eigen_values,eigen_vectors = eigen(LDmatrix_blk)
        # Sort the eigenvalues in decreasing order.
        decreasing_idx = sortperm(eigen_values, rev=true)
        eigen_values      = eigen_values[decreasing_idx]
        eigen_vectors     = eigen_vectors[:,decreasing_idx]
        cumulative_sum    = cumsum(eigen_values)
        # Find the index of the first eigenvalue that is greater than or equal to 99.5% of the total variance.
        stop_index        = findfirst(cumulative_sum .>= 0.995 * sum(eigen_values))
        eigen_value_used     = eigen_values[1:stop_index]
        eigen_vectors_used   = eigen_vectors[:,1:stop_index]
        TransformedX_dict_ranki[blk] = Diagonal(sqrt.(eigen_value_used)) * eigen_vectors_used'
        TransformedY_dict_ranki[blk] = Diagonal(1 ./ sqrt.(eigen_value_used)) * eigen_vectors_used' * bhat_dict_ranki[blk] 
    end
    JLD2.save(mpi_data_path*"rank$rankID.TransformedX_dict.jld2","my_TransformedX_dict",TransformedX_dict_ranki)
    JLD2.save(mpi_data_path*"rank$rankID.TransformedY_dict.jld2", "my_TransformedY_dict",TransformedY_dict_ranki) 
    JLD2.save(mpi_data_path*"rank$rankID.blkSNPsIndex_dict.jld2","my_blkSNPsIndex_dict",blkSNPsIndex_dict_ranki) 
end
