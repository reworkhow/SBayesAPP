using LinearAlgebra, Distributions, Random, SparseArrays
using DelimitedFiles, DataFrames, CSV, JLD2
using Dates
using Statistics,ProgressMeter; 
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--seed"
        help = ""
        arg_type = Int
        default = 123

        "--pleio_percent"
        help = ""
        arg_type = Float64
        default = 0.1

        "--nrank"
        help = ""
        arg_type = Int
        default = 10

        "--sample_size"
        help = ""
        arg_type = Int
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
    end

    return parse_args(s)
end

# Use the parsed arguments
args = parse_commandline()

##############################################################
## Input data
##############################################################
data_folder_path = "/common/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/"
LDinfo_path = "/common/zhao/tianjing/" # /group/qtlchenggrp/jiayi/MTSBayesC/data/UKBioBank/GCTB_ldm/ukbEUR_HM3/
eigen_data_path = "/common/zhao/tianjing/eigen_data/readable_files/"
##############################################################
## Input data
##############################################################


# Assuming `args` is the dictionary returned from the parse_args() function
sample_size = args["sample_size"]
# annotation_size = args["annotation_size"]
pleio_percent = args["pleio_percent"]
if pleio_percent == 0 || pleio_percent == 1
    pleio_percent = Int(pleio_percent) #0, instead of 0.0
end
seed = args["seed"]
h21 = args["h21"] # h2 for trait 1
h22 = args["h22"] # h2 for trait 2
nrank = args["nrank"]

# Optionally, you can print these values to check
println("Sample size: ", sample_size)
# println("Annotation size: ", annotation_size)
println("Pleiotropy percent: ", pleio_percent)
println("Seed: ", seed)
println("Heritability for trait 1 (h21): ", h21)
println("Heritability for trait 2 (h22): ", h22)
println("Number of ranks for analysis: ", nrank)
println("----Check Point A")

# Variables should be assigned or extracted from args before using them in paths
#data_folder_name = "h2_trait1.$h21.h2_trait2.$(h22)_pleioPercent$(pleio_percent)_sampleSize$(sample_size)_annotationSize$(annotation_size)_seed$seed"
data_folder_name = "h2_trait1.$h21.h2_trait2.$(h22)_pleioPercent$(pleio_percent)_sampleSize$(sample_size)_seed$seed"
println("Data folder name: ", data_folder_name)

# Construct the full path by concatenating the folder path and name
gwas_path = "$(data_folder_path)$(data_folder_name)/"  # Make sure to include the final slash for directories
println("----Check Point B")

# Function to perform join, merge order, and sort
function process_bhat_df(snp_ids_blk, bhat, aux_order_df)
    bhat_df = DataFrame(SNP=snp_ids_blk)
    bhat_df = leftjoin(bhat_df, bhat, on=:SNP)
    bhat_df = leftjoin(bhat_df, aux_order_df, on=:SNP)
    sort!(bhat_df, :Order)
    return bhat_df
end
println("----Check Point C")


# read bhat for all SNPs having values 
bhat1 = CSV.read(gwas_path * "standardPhenoVar/b1_complete.txt", DataFrame)
bhat2 = CSV.read(gwas_path * "standardPhenoVar/b2_complete.txt", DataFrame)
# SNPsc1 = vec(readdlm(gwas_path * "SNPc1.txt", ','))
# SNPsc2 = vec(readdlm(gwas_path * "SNPc2.txt", ','))
SNPsc1 = vec(readdlm(gwas_path * "../SNPc1.txt", ','))
SNPsc2 = vec(readdlm(gwas_path * "../SNPc2.txt", ','))
println("----Check Point D")

nBlocks = 49 #591
blkIDs = collect(1:nBlocks)

# Create a dictionary of sub-matrices
blkSNPsIndex_dict = Dict{Int, Vector{Int64}}()
TransformedX_dict = Dict{Int,Matrix{Float64}}()
TransformedY_dict = Dict{Int,Vector{Vector{Float64}}}()
anno_matrix_dict = Dict{Int,Matrix}()
println("----check point 1")


for blk in blkIDs 
    @show blk
    snp_ids_blk = String.(vec(readdlm(LDinfo_path * "SNPsPerBlock/SNPs_block$blk.csv")))
    nSNPs = length(snp_ids_blk)
    blk_indices = collect(1:nSNPs)
    blkSNPsIndex_dict[blk] = blk_indices 

    # Create an auxiliary DataFrame with SNP and its order in snp_ids_blk
    aux_order_df = DataFrame(SNP=snp_ids_blk, Order=1:length(snp_ids_blk))

    bhat_df1 = process_bhat_df(snp_ids_blk, bhat1, aux_order_df)
    bhat_df2 = process_bhat_df(snp_ids_blk, bhat2, aux_order_df)
    bhat_vec = [bhat_df1[!, :bAdj], bhat_df2[!, :bAdj]] # the order of bhat1 and bhat2 should is the same as snp_ids_blk
    

    # transformed w and Q
    # eigen decomposition from GCTB readableFiles
    eigen_values = readdlm(eigen_data_path * "block$blk.lambda.csv")[:, 1]
    eigen_vectors = readdlm(eigen_data_path * "block$blk.U.csv", ',')
    TransformedX_dict[blk] = Diagonal(sqrt.(eigen_values)) * eigen_vectors'
    TransformedY_dict[blk] = [Diagonal(1 ./ sqrt.(eigen_values)) * eigen_vectors' * bhat_vec[b] for b in 1:2]

    # generate Amat dictionary (annot)
    # Make sure the order of annot_df is the same as the order of snp_ids_blk

    anno_matrix_dict_blk = zeros(nSNPs, 2)
    common_SNPsc1 = intersect(snp_ids_blk, SNPsc1)
    common_SNPsc2 = intersect(snp_ids_blk, SNPsc2)

    if length(common_SNPsc1) + length(common_SNPsc2) != length(snp_ids_blk)
        error("annotation file wrong")
    end

    # Find the index of these common elements in snp_ids_blk
    index_common_SNPsc1 = findall(x -> x in common_SNPsc1, snp_ids_blk)
    index_common_SNPsc2 = findall(x -> x in common_SNPsc2, snp_ids_blk)
    anno_matrix_dict_blk[index_common_SNPsc1, 1] .= 1
    anno_matrix_dict_blk[index_common_SNPsc2, 2] .= 1
    anno_matrix_dict[blk] = anno_matrix_dict_blk
end
println("----check point 2")

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
println("----check point 3")

#save data into data folder
mpi_data_path = gwas_path * "nrank$nrank.eigen/bhatXsj/995Eigen/"
mkpath(mpi_data_path)

for i in 1:nrank
    rankID = i - 1 #rank starts from 0, not 1
    writedlm(mpi_data_path * "rank$rankID.blkIDs.txt", blkID_each_rank[i], ',')
end
println("----check point 4")

for i in 1:nrank
    @show i
    rankID = i - 1 #rank starts from 0, not 1
    blkID_ranki = Int.(vec(readdlm(mpi_data_path * "rank$rankID.blkIDs.txt", ',')))  #block ID for this rank

    blkSNPsIndex_dict_ranki = Dict(key => blkSNPsIndex_dict[key] for key in blkID_ranki) #get blkSNPsIndex_dict for this rank
    JLD2.save(mpi_data_path * "rank$rankID.blkSNPsIndex_dict.jld2", "my_blkSNPsIndex_dict", blkSNPsIndex_dict_ranki)  #save

    TransformedX_dict_ranki = Dict(key => TransformedX_dict[key] for key in blkID_ranki)
    JLD2.save(mpi_data_path * "rank$rankID.TransformedX_dict.jld2", "my_TransformedX_dict", TransformedX_dict_ranki)  #save

    TransformedY_dict_ranki = Dict(key => TransformedY_dict[key] for key in blkID_ranki)
    JLD2.save(mpi_data_path * "rank$rankID.TransformedY_dict.jld2", "my_TransformedY_dict", TransformedY_dict_ranki)  #save

    anno_matrix_dict_ranki = Dict(key => anno_matrix_dict[key] for key in blkID_ranki)
    JLD2.save(mpi_data_path * "rank$rankID.anno_matrix_dict.jld2", "my_anno_matrix_dict", anno_matrix_dict_ranki)  #save
end


