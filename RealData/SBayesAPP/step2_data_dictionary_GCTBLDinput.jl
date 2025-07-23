using LinearAlgebra, Distributions, Random, SparseArrays
using DelimitedFiles, DataFrames, CSV, JLD2
using Dates
using Statistics,ProgressMeter; 

data_path = ARGS[1]  # this is the preprocess_file_path in step1 : /group/qtlchenggrp/jiayi/MTSBayesC/data/UKBioBank/RealData/T2D_Comorbidities/T2D_AD/
                     # /common/zhao/jyqqu/MTSBayesCC/data/real_data/T2D_FG/
nrank     = ARGS[2]  # nrank = 20 -> number of ranks to be used in MPI
annot_file = ARGS[3] # annotation file; ../cell_type_annot_human_total_unoverlap.txt
LDinfo_path = ARGS[4] # /group/qtlchenggrp/jiayi/MTSBayesC/data/UKBioBank/GCTB_ldm/ukbEUR_HM3/
                      # /common/zhao/jyqqu/MTSBayesCC/data/eigen_data_whole_genome/
annot_dict_name = ARGS[5] # annotation dictionary name; anno_matrix_cell_type_human_total_unoverlap_dict
nrank     = parse(Int64, nrank)



# Function to perform join, merge order, and sort
function process_bhat_df(snp_ids_blk, bhat, aux_order_df)
    bhat_df = DataFrame(SNP=snp_ids_blk)
    bhat_df = leftjoin(bhat_df, bhat, on=:SNP)
    bhat_df = leftjoin(bhat_df, aux_order_df, on=:SNP)
    sort!(bhat_df, :Order)
    return bhat_df
end


# # read bhat for all SNPs having values 
# bhat1 = CSV.read(data_path * "standardPhenoVar/b1_complete.txt", DataFrame) 
# bhat2 = CSV.read(data_path * "standardPhenoVar/b2_complete.txt", DataFrame)

# annotation df
# first column is SNP, the rest are annotations
annot = CSV.read(data_path * annot_file, DataFrame)


nBlocks = 591
blkIDs = collect(1:nBlocks)

# Create a dictionary of sub-matrices
# blkSNPsIndex_dict = Dict{Int, Vector{Int64}}()
# TransformedX_dict = Dict{Int,Matrix{Float64}}()
# TransformedY_dict = Dict{Int,Vector{Vector{Float64}}}()
# nGWAS_dict = Dict{Int,Vector{Float64}}()
anno_matrix_dict = Dict{Int,Matrix}()


for blk in blkIDs 
    @show blk
    snp_ids_blk = String.(vec(readdlm(LDinfo_path * "SNPsPerBlock/SNPs_block$blk.csv")))
    nSNPs = length(snp_ids_blk)
    # blk_indices = collect(1:nSNPs)
    # blkSNPsIndex_dict[blk] = blk_indices 

    # Create an auxiliary DataFrame with SNP and its order in snp_ids_blk
    aux_order_df = DataFrame(SNP=snp_ids_blk, Order=1:length(snp_ids_blk))

    # bhat_df1 = process_bhat_df(snp_ids_blk, bhat1, aux_order_df)
    # bhat_df2 = process_bhat_df(snp_ids_blk, bhat2, aux_order_df)
    # bhat_vec = [bhat_df1[!, :bAdj], bhat_df2[!, :bAdj]] # the order of bhat1 and bhat2 should is the same as snp_ids_blk
    # nGWAS_dict[blk] = [mean(bhat_df1[!, :N]), mean(bhat_df2[!, :N])]

    # transformed w and Q
    # eigen decomposition from GCTB readableFiles
    # eigen_values = readdlm(LDinfo_path * "readableFiles/block$blk.lambda.csv")[:,1]
    # eigen_vectors = readdlm(LDinfo_path * "readableFiles/block$blk.U.csv", ',')
    # TransformedX_dict[blk] = Diagonal(sqrt.(eigen_values)) * eigen_vectors'
    # TransformedY_dict[blk] = [Diagonal(1 ./ sqrt.(eigen_values)) * eigen_vectors' * bhat_vec[b] for b in 1:2]

    # generate Amat dictionary (annot)
    # Make sure the order of annot_df is the same as the order of snp_ids_blk
    annot_df = DataFrame(SNP=snp_ids_blk)
    annot_df = leftjoin(annot_df, annot, on=:SNP)
    # Replace missing values with 0
    transform!(annot_df, names(annot)[2:end] .=> (x -> coalesce.(x, 0)) .=> names(annot)[2:end])
    annot_df = leftjoin(annot_df, aux_order_df, on=:SNP)
    sort!(annot_df, :Order)
    select!(annot_df, Not(:Order))
    anno_matrix_dict[blk] = Matrix(annot_df[:, 2:end])
end

############################################################
#MPI: split data for different ranks (#rank=#parallel jobs)
############################################################
@show nrank
##############################
# 1. split blkIDs for different ranks
##############################

blkID_each_rank = Vector{Vector{Int}}(undef, nrank)
nPervec = Int(ceil(nBlocks / nrank))
# Split the vector into sub-vectors
for i in 1:nrank
    start_index = (i - 1) * nPervec + 1
    end_index = min(i * nPervec, length(blkIDs))
    blkID_each_rank[i] = blkIDs[start_index:end_index]
end

#save data into data folder
mpi_data_path = data_path * "nrank$nrank.eigen/bhatXsj/995Eigen/"
mkpath(mpi_data_path)

# for i in 1:nrank
#     rankID = i - 1 #rank starts from 0, not 1
#     writedlm(mpi_data_path * "rank$rankID.blkIDs.txt", blkID_each_rank[i], ',')
# end

for i in 1:nrank
    @show i
    rankID = i - 1 #rank starts from 0, not 1
    blkID_ranki = Int.(vec(readdlm(mpi_data_path * "rank$rankID.blkIDs.txt", ',')))  #block ID for this rank

    # blkSNPsIndex_dict_ranki = Dict(key => blkSNPsIndex_dict[key] for key in blkID_ranki) #get blkSNPsIndex_dict for this rank
    # JLD2.save(mpi_data_path * "rank$rankID.blkSNPsIndex_dict.jld2", "my_blkSNPsIndex_dict", blkSNPsIndex_dict_ranki)  #save

    # nGWAS_dict_ranki = Dict(key => nGWAS_dict[key] for key in blkID_ranki)
    # JLD2.save(mpi_data_path * "rank$rankID.nGWAS_dict.jld2", "my_nGWAS_dict", nGWAS_dict_ranki)  #save
    
    # TransformedX_dict_ranki = Dict(key => TransformedX_dict[key] for key in blkID_ranki)
    # JLD2.save(mpi_data_path * "rank$rankID.TransformedX_dict.jld2", "my_TransformedX_dict", TransformedX_dict_ranki)  #save

    # TransformedY_dict_ranki = Dict(key => TransformedY_dict[key] for key in blkID_ranki)
    # JLD2.save(mpi_data_path * "rank$rankID.TransformedY_dict.jld2", "my_TransformedY_dict", TransformedY_dict_ranki)  #save

    anno_matrix_dict_ranki = Dict(key => anno_matrix_dict[key] for key in blkID_ranki)
    JLD2.save(mpi_data_path * "rank$rankID.$annot_dict_name.jld2", "my_anno_matrix_dict", anno_matrix_dict_ranki)  #save
end


