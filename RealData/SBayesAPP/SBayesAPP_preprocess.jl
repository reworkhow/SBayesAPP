using CSV, DataFrames, Statistics, DelimitedFiles
using LinearAlgebra, Distributions, Random, SparseArrays
using JLD2, Dates, ProgressMeter
using ArgParse, Printf

# =================== Command Line Arguments ===================
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--LD_info_path"
            help = "Path to LD info files"
            arg_type = String
        "--out"
            help = "Output path"
            arg_type = String
        "--trait1_file"
            help = "Trait 1 summary file"
            arg_type = String
        "--trait2_file"
            help = "Trait 2 summary file"
            arg_type = String
        "--LDinfo_file"
            help = "LD SNP info file"
            arg_type = String
        "--annot_file"
            help = "Annotation file"
            arg_type = String
        "--annot_dict_name"
            help = "Annotation dict name for JLD2 output"
            arg_type = String
            default = "anno_matrix_dict"
        "--nrank"
            help = "Number of MPI ranks"
            arg_type = Int
        "--nblock"
            help = "Number of blocks to process"
            arg_type = Int
            default = 591
    end
    return parse_args(s)
end

args = parse_commandline()
LD_info_path = args["LD_info_path"]
out_path = args["out"]
trait1_file = args["trait1_file"]
trait2_file = args["trait2_file"]
LDinfo_file = args["LDinfo_file"]
annot_file = args["annot_file"]
annot_dict_name = args["annot_dict_name"]
nrank = args["nrank"]
nBlocks = args["nblock"]

# LD_info_path   = "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/eigen_data_whole_genome/"
# out_path       = "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/h2_trait1.0.5.h2_trait2.0.2_pleioPercent0.5_sampleSize300000_seed1/"
# trait1_file    = "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/h2_trait1.0.5.h2_trait2.0.2_pleioPercent0.5_sampleSize300000_seed1/Trait1.gcta.phen.plink.ci.assoc.linear.ma"
# trait2_file    = "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/h2_trait1.0.5.h2_trait2.0.2_pleioPercent0.5_sampleSize300000_seed1/Trait2.gcta.phen.plink.ci.assoc.linear.ma"
# LDinfo_file    = "snp.info"
# annot_file     = "/mnt/nrdstor/zhao/jyqqu/MTSBayesCC/data/sim_chr1_output_v2/example_analysis/annotation_df.txt"
# annot_dict_name= "anno_matrix_dict"
# nrank = 2
# nBlocks = 49

# =================== Detect Delimiter Function ===================
# This function detects the delimiter of a file based on its first line.
# It returns '\t' for tab, ' ' for space, ',' for comma, or raises an error if the delimiter is not recognized.
function detect_delimiter(file_path::String)
    first_line = strip(readline(file_path))

    if occursin(r"\t", first_line)
        return ('\t', false)
    elseif occursin(r"  +", first_line)  # two or more spaces
        return (' ', true)
    elseif occursin(r"^[^\t ]+( [^\t ]+)+$", first_line)  # single-space separation
        return (' ', false)
    elseif occursin(',', first_line)
        return (',', false)
    else
        error("Unknown delimiter format.")
    end
end

# =================== Adjust Effect Size ===================
function adjust_effects!(bhat::DataFrame, blockinfo::DataFrame)
    bhat.freq = ifelse.(bhat.A1 .== blockinfo.A2, 1 .- bhat.freq, bhat.freq)
    bhat.b    = ifelse.(bhat.A1 .== blockinfo.A2, -bhat.b, bhat.b)
    bhat.A1   = blockinfo.A1
    bhat.A2   = blockinfo.A2

    N = bhat.N
    se = bhat.se
    b = bhat.b
    
    sj = sqrt.(1.0 ./ (N .* se.^2 .+ b.^2))
    bhat.bAdj = b .* sj
    bhat.seAdj = se .* sj
    bhat.D = 2 .* bhat.freq .* (1 .- bhat.freq) .* N
    bhat.varps = N .* bhat.seAdj.^2 .+ bhat.bAdj.^2

    return bhat
end

# =================== Order bhat based LD matrix order ===================
function process_bhat_df(snp_ids_blk, bhat, aux_order_df)
    bhat_df = DataFrame(SNP=snp_ids_blk)
    bhat_df = leftjoin(bhat_df, bhat, on=:SNP)
    bhat_df = leftjoin(bhat_df, aux_order_df, on=:SNP)
    sort!(bhat_df, :Order)
    return bhat_df
end

# =================== Load Data ===================
delim, ignore_repeated = detect_delimiter(trait1_file)
bhat1 = CSV.read(trait1_file, DataFrame; delim=delim, ignorerepeated=ignore_repeated)
bhat2 = CSV.read(trait2_file, DataFrame; delim=delim, ignorerepeated=ignore_repeated)
delim, ignore_repeated = detect_delimiter(joinpath(LD_info_path, LDinfo_file))
LDinfo = CSV.read(joinpath(LD_info_path, LDinfo_file), DataFrame; delim=delim, ignorerepeated=ignore_repeated)
delim, ignore_repeated = detect_delimiter(annot_file)
annot = CSV.read(annot_file, DataFrame; delim=delim, ignorerepeated=ignore_repeated)

mkpath(joinpath(LD_info_path, "SNPsPerBlock_test"))
mkpath(joinpath(out_path, "standardPhenoVar"))

# =================== Loop Over LD Blocks and Standardize ===================
println("Standardizing SNP marginal effects...")
for i in 1:nBlocks
    blockinfo = filter(:Block => x -> x == i, LDinfo)
    snps_ld_order = blockinfo.ID

    bhat1i = bhat1[in.(bhat1.SNP, Ref(blockinfo.ID)), :]
    bhat2i = bhat2[in.(bhat2.SNP, Ref(blockinfo.ID)), :]

    # Match bhat1i / bhat2i order to blockinfo.ID
    bhat1i = leftjoin(blockinfo[:, [:ID, :Index]], bhat1i, on = [:ID => :SNP])
    rename!(bhat1i, :ID => :SNP)   # Rename ID back to SNP
    sort!(bhat1i, :Index)
    bhat2i = leftjoin(blockinfo[:, [:ID, :Index]], bhat2i, on = [:ID => :SNP])
    rename!(bhat2i, :ID => :SNP)   # Rename ID back to SNP
    sort!(bhat2i, :Index)

    if !(all(blockinfo.ID .== bhat1i.SNP) && all(blockinfo.ID .== bhat2i.SNP))
        error("Mismatch in SNP order for block $i")
    end

    adjust_effects!(bhat1i, blockinfo)
    adjust_effects!(bhat2i, blockinfo)

    # Add empty columns if they do not exist
    for col in ["bAdj", "varps"]
        for df in [bhat1, bhat2]
            if !(col in names(df))
                df[!, Symbol(col)] = Vector{Union{Missing, Float64}}(missing, nrow(df))
            end
        end
    end

    for (bhat, trait) in [(bhat1i, bhat1), (bhat2i, bhat2)]
        for row in eachrow(bhat)
            idx = findfirst(==(row.SNP), trait.SNP)
            if isnothing(idx); continue; end
            trait[idx, [:bAdj, :varps, :A1, :A2, :b, :freq]] = row[[:bAdj, :varps, :A1, :A2, :b, :freq]]
        end
    end

    writedlm(joinpath(LD_info_path, "SNPsPerBlock_test", "SNPs_block$(i).csv"), snps_ld_order)
end
# ========== Save standardized marginal effects ==========
CSV.write(joinpath(out_path, "standardPhenoVar", "b1_complete.txt"), bhat1)
CSV.write(joinpath(out_path, "standardPhenoVar", "b2_complete.txt"), bhat2)

# =================== Process MPI Data for Blocks ====================
println("Processing blocks for MPI...")
blkIDs = collect(1:nBlocks)
# Create dictionaries to save sub-matrices
blkSNPsIndex_dict = Dict{Int, Vector{Int64}}()
TransformedX_dict = Dict{Int, Matrix}()
TransformedY_dict = Dict{Int, Vector{Vector{Float64}}}()
nGWAS_dict = Dict{Int, Vector{Float64}}()
anno_matrix_dict = Dict{Int, Matrix}()

for blk in blkIDs
    snp_ids_blk = String.(vec(readdlm(joinpath(LD_info_path, "SNPsPerBlock_test", "SNPs_block$blk.csv"))))
    blkSNPsIndex_dict[blk] = collect(1:length(snp_ids_blk))
    aux_order_df = DataFrame(SNP=snp_ids_blk, Order=1:length(snp_ids_blk))

    bhat_df1 = process_bhat_df(snp_ids_blk, bhat1, aux_order_df)
    bhat_df2 = process_bhat_df(snp_ids_blk, bhat2, aux_order_df)
    bhat_vec = [bhat_df1.bAdj, bhat_df2.bAdj]
    nGWAS_dict[blk] = [mean(bhat_df1.N), mean(bhat_df2.N)]

    λ = readdlm(joinpath(LD_info_path, "readableFiles", "block$blk.lambda.csv"))[:, 1]
    U = readdlm(joinpath(LD_info_path, "readableFiles", "block$blk.U.csv"), ',')

    TransformedX_dict[blk] = Diagonal(sqrt.(λ)) * U'
    TransformedY_dict[blk] = [Diagonal(1 ./ sqrt.(λ)) * U' * bhat_vec[b] for b in 1:2]

    # Generate annotation matrix
    annot_df = leftjoin(DataFrame(SNP=snp_ids_blk), annot, on=:SNP)
    # Replace missing values with 0
    transform!(annot_df, names(annot)[2:end] .=> (x -> coalesce.(x, 0)) .=> names(annot)[2:end])
    annot_df = leftjoin(annot_df, aux_order_df, on=:SNP)
    sort!(annot_df, :Order)
    select!(annot_df, Not(:Order))
    anno_matrix_dict[blk] = Matrix(annot_df[:, 2:end])
end

println("Saving MPI data in each rank...")
# Split & Save for Each MPI Rank
blkID_each_rank = Vector{Vector{Int}}(undef, nrank)
nPervec = Int(ceil(nBlocks / nrank))
# Split the vector into sub-vectors
for i in 1:nrank
    start_index = (i - 1) * nPervec + 1
    end_index = min(i * nPervec, length(blkIDs))
    blkID_each_rank[i] = blkIDs[start_index:end_index]
end
#save data into out folder
mpi_path = joinpath(out_path, "nrank$nrank.eigen", "bhatXsj", "995Eigen")
mkpath(mpi_path)

for i in 1:nrank
    rankID = i - 1
    blk_set = blkID_each_rank[i]
    writedlm(joinpath(mpi_path, "rank$rankID.blkIDs.txt"), blk_set, ',')

    JLD2.save(joinpath(mpi_path, "rank$rankID.blkSNPsIndex_dict.jld2"), "my_blkSNPsIndex_dict", Dict(k => blkSNPsIndex_dict[k] for k in blk_set))
    JLD2.save(joinpath(mpi_path, "rank$rankID.nGWAS_dict.jld2"), "my_nGWAS_dict", Dict(k => nGWAS_dict[k] for k in blk_set))
    JLD2.save(joinpath(mpi_path, "rank$rankID.TransformedX_dict.jld2"), "my_TransformedX_dict", Dict(k => TransformedX_dict[k] for k in blk_set))
    JLD2.save(joinpath(mpi_path, "rank$rankID.TransformedY_dict.jld2"), "my_TransformedY_dict", Dict(k => TransformedY_dict[k] for k in blk_set))
    JLD2.save(joinpath(mpi_path, "rank$rankID.$annot_dict_name.jld2"), "my_anno_matrix_dict", Dict(k => anno_matrix_dict[k] for k in blk_set))
end
println("Preprocessing completed successfully!")
# ========== End of Script ==========

