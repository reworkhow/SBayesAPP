using CSV
using DataFrames
using DelimitedFiles: readdlm, writedlm
using JLD2
using LinearAlgebra: Diagonal
using Statistics: mean

mutable struct NonMPIBlockData{TX,TY,TS,TN,TA}
    block_id::Int
    transformed_x::TX
    transformed_y::TY
    snp_indices::TS
    n_gwas::TN
    anno_matrix::TA
    xpx::Union{Nothing,Vector{Vector{Float64}}}
    x_arrays::Union{Nothing,Vector{Matrix{Float64}}}
    annotation_mask::Union{Nothing,AbstractMatrix{Bool}}
end

struct NonMPIBlockCollection{TB}
    blocks::Vector{TB}
    nblk::Int
    nsnp::Int
end

const REQUIRED_GCTB_MA_COLUMNS = ("SNP", "A1", "A2", "freq", "b", "se", "N")
const REQUIRED_LDINFO_COLUMNS = ("Block", "ID", "Index", "A1", "A2")

function validate_gctb_ma_columns(df::DataFrame, data_name::AbstractString)
    missing_columns = [column for column in REQUIRED_GCTB_MA_COLUMNS if !(column in names(df))]
    isempty(missing_columns) || error(
        "$data_name must contain GCTB .ma columns $(join(REQUIRED_GCTB_MA_COLUMNS, ", ")). Missing: $(join(missing_columns, ", "))",
    )
    return nothing
end

function validate_ldinfo_columns(df::DataFrame, data_name::AbstractString)
    missing_columns = [column for column in REQUIRED_LDINFO_COLUMNS if !(column in names(df))]
    isempty(missing_columns) || error(
        "$data_name must contain LD info columns $(join(REQUIRED_LDINFO_COLUMNS, ", ")). Missing: $(join(missing_columns, ", "))",
    )
    return nothing
end

function load_jld2_entry(path::AbstractString, candidate_keys)
    data = JLD2.load(path)
    for key in candidate_keys
        haskey(data, key) && return data[key]
    end
    error("None of keys $(collect(candidate_keys)) found in $path. Available keys: $(collect(keys(data)))")
end

function load_annotation_metadata(data_path::AbstractString, annot_file::AbstractString; nCon::Int=0)
    annot = CSV.read(data_path * annot_file, DataFrame)
    annotation_name = String.(names(annot)[2:end])
    n_annotation = length(annotation_name)
    0 <= nCon <= n_annotation || error("nCon must be between 0 and $n_annotation for $(data_path * annot_file)")

    annotation_values = Matrix(annot[!, 2:end])
    n_loci_annot = Int.(vec(sum(annotation_values .!= 0, dims=1)))
    n_cat = n_annotation - nCon
    annotation_type = vcat(fill("continue", nCon), fill("category", n_cat))

    return (
        annotationName=annotation_name,
        nLoci_annot=n_loci_annot,
        nCon=nCon,
        nCat=n_cat,
        annotationType=annotation_type,
    )
end

function detect_input_delimiter(path::AbstractString)
    return open(path, "r") do io
        eof(io) && error("No header found in $path")
        first_line = rstrip(readline(io), ['\r', '\n'])
        if length(split(first_line, '\t')) > 1
            return ('\t', false)
        elseif length(split(first_line, ',')) > 1
            return (',', false)
        elseif length(split(first_line)) > 1
            return (' ', true)
        end
        error("Unknown delimiter format in $path")
    end
end

function read_input_table(path::AbstractString)
    delim, ignore_repeated = detect_input_delimiter(path)
    return CSV.read(path, DataFrame; delim=delim, ignorerepeated=ignore_repeated)
end

function resolve_input_path(data_path::AbstractString, path::AbstractString)
    candidate = String(path)
    isabspath(candidate) && return candidate
    isfile(candidate) && return abspath(candidate)
    return joinpath(data_path, candidate)
end

default_annotation_dict_filename(annot_path::AbstractString) = "anno_matrix_" * splitext(basename(annot_path))[1] * ".jld2"

normalize_dir_path(path::AbstractString) = endswith(path, "/") ? String(path) : string(path, "/")

function align_trait_to_ldinfo(ldinfo::DataFrame, trait_df::DataFrame, trait_label::AbstractString)
    aligned_trait = leftjoin(
        select(ldinfo, :Block, :ID, :Index, :A1 => :LD_A1, :A2 => :LD_A2),
        trait_df,
        on=[:ID => :SNP],
    )
    sort!(aligned_trait, [:Block, :Index])

    missing_mask = ismissing.(aligned_trait[!, :A1])
    if any(missing_mask)
        missing_snps = String.(aligned_trait[missing_mask, :ID])
        n_missing_snps = length(missing_snps)
        error(
            "$trait_label is missing $n_missing_snps SNPs (based on A1) from LD info. First missing SNPs: $(join(missing_snps[1:min(n_missing_snps, 5)], ", "))",
        )
    end

    return aligned_trait
end

function standardize_summary_block(
    block_trait::DataFrame;
    block_id::Integer,
    trait_label::AbstractString,
)
    trait_a1 = String.(block_trait[!, :A1])
    trait_a2 = String.(block_trait[!, :A2])
    block_a1 = String.(block_trait[!, :LD_A1])
    block_a2 = String.(block_trait[!, :LD_A2])

    allele_match = ((trait_a1 .== block_a1) .& (trait_a2 .== block_a2)) .| ((trait_a1 .== block_a2) .& (trait_a2 .== block_a1))
    if !all(allele_match)
        mismatch_rows = findall(.!allele_match)
        error(
            "$trait_label has allele mismatches against LD info in block $block_id. Mismatch row indices: $(join(mismatch_rows, ", "))",
        )
    end

    freq = Float64.(block_trait[!, :freq])
    effect = Float64.(block_trait[!, :b])
    se = Float64.(block_trait[!, :se])
    n = Float64.(block_trait[!, :N])
    flip = trait_a1 .== block_a2

    adjusted_freq = ifelse.(flip, 1.0 .- freq, freq)
    adjusted_b = ifelse.(flip, -effect, effect)
    sj = sqrt.(1.0 ./ (n .* se.^2 .+ adjusted_b.^2))
    adjusted_block = copy(block_trait)
    adjusted_block[!, :A1] = block_a1
    adjusted_block[!, :A2] = block_a2
    adjusted_block[!, :freq] = adjusted_freq
    adjusted_block[!, :b] = adjusted_b
    adjusted_block[!, :bAdj] = adjusted_b .* sj
    adjusted_block[!, :seAdj] = se .* sj
    adjusted_block[!, :D] = 2.0 .* adjusted_freq .* (1.0 .- adjusted_freq) .* n
    adjusted_block[!, :varps] = n .* adjusted_block[!, :seAdj].^2 .+ adjusted_block[!, :bAdj].^2

    return adjusted_block
end

function write_jld2_dict(path::AbstractString, key::AbstractString, value)
    jldopen(path, "w") do file
        write(file, key, value)
    end
    return path
end

function build_nonmpi_input_dicts(
    output_path::AbstractString,
    ld_info_path::AbstractString,
    trait1_file::AbstractString,
    trait2_file::AbstractString,
    ldinfo_file::AbstractString,
    annot_file::AbstractString;
    annot_dict_name::AbstractString="anno_matrix_dict",
    readable_files_dir::AbstractString="readableFiles",
    nblocks::Union{Nothing,Integer}=nothing,
)
    output_dir = normalize_dir_path(isabspath(output_path) ? output_path : abspath(output_path))
    ld_root = normalize_dir_path(isabspath(ld_info_path) ? ld_info_path : abspath(ld_info_path))
    trait1_path = resolve_input_path(pwd(), trait1_file)
    trait2_path = resolve_input_path(pwd(), trait2_file)
    ldinfo_path = resolve_input_path(ld_root, ldinfo_file)
    annot_path = resolve_input_path(pwd(), annot_file)
    readable_dir = joinpath(ld_root, readable_files_dir)

    mkpath(output_dir)
    isdir(readable_dir) || error("Readable LD directory does not exist: $readable_dir")

    trait1 = read_input_table(trait1_path)
    trait2 = read_input_table(trait2_path)
    ldinfo = read_input_table(ldinfo_path)

    validate_gctb_ma_columns(trait1, "Trait 1 summary")
    validate_gctb_ma_columns(trait2, "Trait 2 summary")
    validate_ldinfo_columns(ldinfo, "LD info")
    sort!(ldinfo, [:Block, :Index])

    trait1_aligned = align_trait_to_ldinfo(ldinfo, trait1, "Trait 1 summary")
    trait2_aligned = align_trait_to_ldinfo(ldinfo, trait2, "Trait 2 summary")

    all_block_ids = sort(unique(Int.(ldinfo[!, :Block])))
    block_ids = isnothing(nblocks) ? all_block_ids : first(all_block_ids, min(length(all_block_ids), Int(nblocks)))

    transformed_x_dict = Dict{Int,Matrix{Float64}}()
    transformed_y_dict = Dict{Int,Vector{Vector{Float64}}}()
    blk_snps_index_dict = Dict{Int,Vector{Int64}}()
    nGWAS_dict = Dict{Int,Vector{Float64}}()
    adjusted_trait1_blocks = DataFrame[]
    adjusted_trait2_blocks = DataFrame[]

    for block_id in block_ids
        blockinfo = ldinfo[Int.(ldinfo[!, :Block]) .== block_id, :]
        trait1_block_df = trait1_aligned[Int.(trait1_aligned[!, :Block]) .== block_id, :]
        trait2_block_df = trait2_aligned[Int.(trait2_aligned[!, :Block]) .== block_id, :]

        block_snps = String.(blockinfo[!, :ID])
        all(block_snps .== String.(trait1_block_df[!, :ID])) || error("Trait 1 summary SNP order mismatch for block $block_id")
        all(block_snps .== String.(trait2_block_df[!, :ID])) || error("Trait 2 summary SNP order mismatch for block $block_id")

        trait1_block = standardize_summary_block(
            trait1_block_df;
            block_id=block_id,
            trait_label="Trait 1 summary",
        )
        trait2_block = standardize_summary_block(
            trait2_block_df;
            block_id=block_id,
            trait_label="Trait 2 summary",
        )
        push!(adjusted_trait1_blocks, trait1_block)
        push!(adjusted_trait2_blocks, trait2_block)

        lambda_path = joinpath(readable_dir, "block$(block_id).lambda.csv")
        u_path = joinpath(readable_dir, "block$(block_id).U.csv")
        isfile(lambda_path) || error("Missing eigenvalue file for block $block_id: $lambda_path")
        isfile(u_path) || error("Missing eigenvector file for block $block_id: $u_path")

        lambda_raw = readdlm(lambda_path)
        lambda = Float64.(lambda_raw[:, 1])
        u_raw = readdlm(u_path, ',')
        u_matrix = Matrix{Float64}(u_raw)

        length(lambda) == size(u_matrix, 2) || error("Block $block_id has $(length(lambda)) eigenvalues but $(size(u_matrix, 2)) eigenvector columns.")
        size(u_matrix, 1) == length(block_snps) || error("Block $block_id has $(size(u_matrix, 1)) eigenvector rows but $(length(block_snps)) SNPs in LD info.")

        sqrt_lambda_u = Diagonal(sqrt.(lambda)) * transpose(u_matrix)
        inv_sqrt_lambda_u = Diagonal(1.0 ./ sqrt.(lambda)) * transpose(u_matrix)

        transformed_x_dict[block_id] = Matrix{Float64}(sqrt_lambda_u)
        transformed_y_dict[block_id] = [
            Vector{Float64}(inv_sqrt_lambda_u * trait1_block.bAdj),
            Vector{Float64}(inv_sqrt_lambda_u * trait2_block.bAdj),
        ]
        blk_snps_index_dict[block_id] = collect(Int64(1):Int64(length(block_snps)))
        nGWAS_dict[block_id] = [mean(Float64.(trait1_block[!, :N])), mean(Float64.(trait2_block[!, :N]))]
    end

    adjusted_output_dir = joinpath(output_dir, "standardPhenoVar")
    mkpath(adjusted_output_dir)
    trait1_adjusted = vcat(adjusted_trait1_blocks...)
    trait2_adjusted = vcat(adjusted_trait2_blocks...)
    trait1_adjusted_file = joinpath(adjusted_output_dir, "trait1_complete.csv")
    trait2_adjusted_file = joinpath(adjusted_output_dir, "trait2_complete.csv")
    CSV.write(trait1_adjusted_file, trait1_adjusted)
    CSV.write(trait2_adjusted_file, trait2_adjusted)

    jld2_output_dir = joinpath(output_dir, "combined_dict")
    mkpath(jld2_output_dir)
    writedlm(joinpath(jld2_output_dir, "blkIDs.txt"), block_ids, ',')
    write_jld2_dict(joinpath(jld2_output_dir, "TransformedX_dict.jld2"), "my_TransformedX_dict", transformed_x_dict)
    write_jld2_dict(joinpath(jld2_output_dir, "TransformedY_dict.jld2"), "my_TransformedY_dict", transformed_y_dict)
    write_jld2_dict(joinpath(jld2_output_dir, "blkSNPsIndex_dict.jld2"), "my_blkSNPsIndex_dict", blk_snps_index_dict)
    write_jld2_dict(joinpath(jld2_output_dir, "nGWAS_dict.jld2"), "my_nGWAS_dict", nGWAS_dict)

    annot_result = build_annotation_dict(
        jld2_output_dir,
        annot_path,
        ldinfo_path;
        output_file_name="$(annot_dict_name).jld2",
    )

    return (
        output_path=output_dir,
        trait1_file=trait1_path,
        trait2_file=trait2_path,
        trait1_adjusted=trait1_adjusted,
        trait2_adjusted=trait2_adjusted,
        trait1_adjusted_file=trait1_adjusted_file,
        trait2_adjusted_file=trait2_adjusted_file,
        ldinfo_file=ldinfo_path,
        annotation_file=annot_path,
        annot_dict=annot_dict_name,
        block_ids=block_ids,
        transformed_x_file=joinpath(jld2_output_dir, "TransformedX_dict.jld2"),
        transformed_y_file=joinpath(jld2_output_dir, "TransformedY_dict.jld2"),
        blk_snps_index_file=joinpath(jld2_output_dir, "blkSNPsIndex_dict.jld2"),
        nGWAS_file=joinpath(jld2_output_dir, "nGWAS_dict.jld2"),
        annotation_dict_file=annot_result.output_file,
    )
end

function build_annotation_dict(
    output_path::AbstractString,
    annot_path::AbstractString,
    ldinfo_path::AbstractString;
    output_file_name::Union{Nothing,AbstractString}=nothing,
)
    output_dir = isabspath(output_path) ? String(output_path) : abspath(output_path)
    blkid_path = joinpath(output_dir, "blkIDs.txt")

    output_file_path = if output_file_name !== nothing
        output_candidate = String(output_file_name)
        isabspath(output_candidate) ? output_candidate : joinpath(output_dir, output_candidate)
    else
        joinpath(output_dir, default_annotation_dict_filename(annot_path))
    end

    annot = read_input_table(annot_path)
    "SNP" in names(annot) || error("Annotation file must contain a SNP column: $annot_path")
    annotation_columns = names(annot)[2:end]
    isempty(annotation_columns) && error("Annotation file must contain at least one annotation column: $annot_path")

    ldinfo = read_input_table(ldinfo_path)
    validate_ldinfo_columns(ldinfo, "LD info")
    sort!(ldinfo, [:Block, :Index])

    block_ids = if isfile(blkid_path)
        sort(Int.(vec(readdlm(blkid_path, ','))))
    else
        sort(unique(Int.(ldinfo[!, :Block])))
    end

    anno_matrix_dict = Dict{Int,Matrix{Float64}}()
    for block_id in block_ids
        blockinfo = ldinfo[Int.(ldinfo[!, :Block]) .== block_id, :]

        block_snps = DataFrame(_row=1:nrow(blockinfo), SNP=String.(blockinfo[!, :ID]))
        annot_df = leftjoin(block_snps, annot, on=:SNP)
        sort!(annot_df, :_row)
        for column in annotation_columns
            annot_df[!, column] = Float64.(coalesce.(annot_df[!, column], 0))
        end

        anno_matrix_dict[block_id] = Matrix(annot_df[:, annotation_columns])
    end

    output_dict_name = "my_anno_matrix_dict"
    jldopen(output_file_path, "w") do file
        write(file, output_dict_name, anno_matrix_dict)
    end

    return (
        output_file=output_file_path,
        output_dict_name=output_dict_name,
        nblocks=length(anno_matrix_dict),
    )
end

function load_nonmpi_block_data(data_path::AbstractString, annot_dict::AbstractString)
    base_dir = normalize_dir_path(isabspath(data_path) ? data_path : abspath(data_path))
    combined_dict_dir = normalize_dir_path(joinpath(base_dir, "combined_dict"))
    dict_dir = isdir(combined_dict_dir) ? combined_dict_dir : base_dir

    transformed_x_dict = load_jld2_entry(dict_dir * "TransformedX_dict.jld2", ("my_TransformedX_dict", "TransformedX_dict"))
    transformed_y_dict = load_jld2_entry(dict_dir * "TransformedY_dict.jld2", ("my_TransformedY_dict", "TransformedY_dict"))
    blkSNPsIndex_dict = load_jld2_entry(dict_dir * "blkSNPsIndex_dict.jld2", ("my_blkSNPsIndex_dict", "blkSNPsIndex_dict"))
    blkID = Int.(vec(readdlm(dict_dir * "blkIDs.txt", ',')))
    nGWAS_dict = load_jld2_entry(dict_dir * "nGWAS_dict.jld2", ("my_nGWAS_dict", "nGWAS_dict"))
    anno_matrix_dict = load_jld2_entry(base_dir * "$annot_dict.jld2", ("my_anno_matrix_dict", annot_dict))
    sort!(blkID)

    blocks = NonMPIBlockData[]
    nsnp = 0
    offset = 0
    for blk in blkID
        raw_indices = Int.(vec(blkSNPsIndex_dict[blk]))
        local_count = length(raw_indices)
        snp_indices = raw_indices .+ offset
        push!(
            blocks,
            NonMPIBlockData(
                blk,
                transformed_x_dict[blk],
                transformed_y_dict[blk],
                snp_indices,
                nGWAS_dict[blk],
                anno_matrix_dict[blk],
                nothing,
                nothing,
                nothing,
            ),
        )
        nsnp += local_count
        offset += local_count
    end

    return NonMPIBlockCollection(blocks, length(blocks), nsnp)
end

function load_effect_state(effect_starting_path, delta_starting_path, my_rank, nTraits)
    betaArray = [
        vec(readdlm(effect_starting_path * "last_mcmc_betaArray$(trait).rank$my_rank.txt"))
        for trait in 1:nTraits
    ]
    deltaArray = [
        vec(readdlm(delta_starting_path * "last_sample_delta$(trait)_rank$my_rank.txt"))
        for trait in 1:nTraits
    ]
    alphaArray = [deltaArray[trait] .* betaArray[trait] for trait in 1:nTraits]
    return betaArray, alphaArray, deltaArray
end