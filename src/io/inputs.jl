using CSV
using DataFrames
using DelimitedFiles: readdlm
using JLD2

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

function build_annotation_dict(
    data_path::AbstractString,
    annot_file::AbstractString,
    ldinfo_file::AbstractString;
    output_file_name::Union{Nothing,AbstractString}=nothing,
)
    annot_path = resolve_input_path(data_path, annot_file)
    ldinfo_path = resolve_input_path(data_path, ldinfo_file)
    blkid_path = joinpath(data_path, "blkIDs.txt")

    output_path = if output_file_name !== nothing
        output_candidate = String(output_file_name)
        isabspath(output_candidate) ? output_candidate : joinpath(data_path, output_candidate)
    else
        joinpath(data_path, default_annotation_dict_filename(annot_path))
    end

    annot = read_input_table(annot_path)
    "SNP" in names(annot) || error("Annotation file must contain a SNP column: $annot_path")
    annotation_columns = names(annot)[2:end]
    isempty(annotation_columns) && error("Annotation file must contain at least one annotation column: $annot_path")

    ldinfo = read_input_table(ldinfo_path)
    "Block" in names(ldinfo) || error("LD info file must contain a Block column: $ldinfo_path")
    "ID" in names(ldinfo) || error("LD info file must contain an ID column: $ldinfo_path")

    block_ids = if isfile(blkid_path)
        sort(Int.(vec(readdlm(blkid_path, ','))))
    else
        sort(unique(Int.(ldinfo[!, :Block])))
    end

    anno_matrix_dict = Dict{Int,Matrix{Float64}}()
    for block_id in block_ids
        blockinfo = ldinfo[Int.(ldinfo[!, :Block]) .== block_id, :]
        "Index" in names(blockinfo) && sort!(blockinfo, :Index)

        block_snps = DataFrame(_row=1:nrow(blockinfo), SNP=String.(blockinfo[!, :ID]))
        annot_df = leftjoin(block_snps, annot, on=:SNP)
        sort!(annot_df, :_row)
        for column in annotation_columns
            annot_df[!, column] = Float64.(coalesce.(annot_df[!, column], 0))
        end

        anno_matrix_dict[block_id] = Matrix(annot_df[:, annotation_columns])
    end

    output_dict_name = "my_anno_matrix_dict"
    jldopen(output_path, "w") do file
        write(file, output_dict_name, anno_matrix_dict)
    end

    return (
        output_file=output_path,
        output_dict_name=output_dict_name,
        nblocks=length(anno_matrix_dict),
    )
end

function load_nonmpi_block_data(data_path::AbstractString, annot_dict::AbstractString)
    transformed_x_dict = load_jld2_entry(data_path * "TransformedX_dict.jld2", ("my_TransformedX_dict", "TransformedX_dict"))
    transformed_y_dict = load_jld2_entry(data_path * "TransformedY_dict.jld2", ("my_TransformedY_dict", "TransformedY_dict"))
    blkSNPsIndex_dict = load_jld2_entry(data_path * "blkSNPsIndex_dict.jld2", ("my_blkSNPsIndex_dict", "blkSNPsIndex_dict"))
    blkID = Int.(vec(readdlm(data_path * "blkIDs.txt", ',')))
    nGWAS_dict = load_jld2_entry(data_path * "nGWAS_dict.jld2", ("my_nGWAS_dict", "nGWAS_dict"))
    anno_matrix_dict = load_jld2_entry(data_path * "$annot_dict.jld2", ("my_anno_matrix_dict", annot_dict))
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