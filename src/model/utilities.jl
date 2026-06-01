using Distributions: InverseWishart
using LinearAlgebra: Symmetric, diag

function parse_cli_args(args::Vector{String}; usage::Union{Nothing,Function}=nothing, allow_empty::Bool=true)
    (args == ["--help"] || args == ["-h"]) && begin
        usage === nothing || usage()
        exit(0)
    end

    isempty(args) && !allow_empty && begin
        usage === nothing || usage()
        exit(1)
    end

    length(args) % 2 == 0 || error("Expected --key value pairs")

    options = Dict{String,String}()
    for index in 1:2:length(args)
        key = args[index]
        value = args[index + 1]
        startswith(key, "--") || error("Expected option starting with --, got: $key")
        options[key[3:end]] = value
    end
    return options
end

option_value(options::AbstractDict, key::AbstractString; default=nothing, required::Bool=false) =
    haskey(options, key) ? options[key] : required ? error("Missing required option --$key") : default

function parse_bool(text::AbstractString)
    lowered = lowercase(strip(text))
    lowered == "true" && return true
    lowered == "false" && return false
    error("Expected true or false, got: $text")
end

function compute_correlation(cov_matrix)
    std_devs = sqrt.(diag(cov_matrix))
    return cov_matrix[1, 2] / (std_devs[1] * std_devs[2])
end

function sample_variance_sumstats(ycorr_array, nobs, df, scale)
    ntraits = length(ycorr_array)
    SSE = zeros(ntraits, ntraits)
    for traiti = 1:ntraits
        ycorri = ycorr_array[traiti]
        for traitj = traiti:ntraits
            ycorrj = ycorr_array[traitj]
            SSE[traiti, traitj] = dot(ycorri, ycorrj)
            SSE[traitj, traiti] = SSE[traiti, traitj]
        end
    end
    return rand(InverseWishart(df + nobs, convert(Array, Symmetric(scale + SSE))))
end

function block_thread_weight(block)
    annotation_mask = getproperty(block, :annotation_mask)
    annotation_mask === nothing && error("block_thread_weight requires prepared block.annotation_mask")
    return max(count(annotation_mask), 1)
end

function build_block_ranges(blocks, ntasks::Int)
    nblocks = length(blocks)
    nblocks == 0 && return UnitRange{Int}[]

    ntasks = max(min(ntasks, nblocks), 1)
    block_weights = [block_thread_weight(block) for block in blocks]
    total_weight = sum(block_weights)
    target_weight = max(cld(total_weight, ntasks), 1)

    ranges = UnitRange{Int}[]
    start_index = 1
    accumulated_weight = 0

    for block_index in 1:nblocks
        accumulated_weight += block_weights[block_index]
        blocks_remaining = nblocks - block_index
        tasks_remaining = ntasks - length(ranges) - 1
        should_split = tasks_remaining > 0 && (
            accumulated_weight >= target_weight ||
            blocks_remaining == tasks_remaining
        )
        if should_split
            push!(ranges, start_index:block_index)
            start_index = block_index + 1
            accumulated_weight = 0
        end
    end

    start_index <= nblocks && push!(ranges, start_index:nblocks)
    return ranges
end

function merge_nloci_counts!(dest, src)
    for cat in eachindex(dest, src)
        dest[cat] .+= src[cat]
    end
    return nothing
end

function two_trait_state_index(delta1::Float64, delta2::Float64)
    if delta1 == 0.0
        return delta2 == 0.0 ? 1 : 4
    end
    return delta2 == 1.0 ? 2 : 3
end

pi_key_order() = ((0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0))

pi_key(key::AbstractVector{<:Real}) = (Float64(key[1]), Float64(key[2]))
pi_key(key::Tuple{Vararg{<:Real,2}}) = (Float64(key[1]), Float64(key[2]))

function format_pi_key(key::Union{AbstractVector{<:Real},Tuple{Vararg{<:Real,2}}})
    return "[" * join(string.(pi_key(key)), ",") * "]"
end

function parse_pi_key(key_text::AbstractString)::NTuple{2,Float64}
    stripped_key = strip(key_text)
    startswith(stripped_key, "[") || error("Invalid Pi key format: $key_text")
    endswith(stripped_key, "]") || error("Invalid Pi key format: $key_text")

    inner_key = stripped_key[2:end-1]
    parts = split(inner_key, ',')
    length(parts) == 2 || error("Invalid Pi key format: $key_text")
    return (parse(Float64, strip(parts[1])), parse(Float64, strip(parts[2])))
end

function split_pi_line(line::AbstractString)
    if occursin('\t', line)
        return split(line, '\t'; limit=2)
    end

    if occursin("],", line)
        key_text, value_text = split(line, "],"; limit=2)
        return [key_text * "]", value_text]
    end

    error("Unsupported Pi file line format: $line")
end

function write_pi_dict(io::IO, dict::AbstractDict)
    for key in pi_key_order()
        haskey(dict, key) || error("Missing Pi key $(format_pi_key(key)) in dictionary output")
        println(io, format_pi_key(key), '\t', dict[key])
    end
    return nothing
end

function write_pi_dict(output_file::AbstractString, dict::AbstractDict)
    open(output_file, "w") do io
        write_pi_dict(io, dict)
    end
    return nothing
end

function read_to_dict(input_file::String)::Dict{NTuple{2,Float64},Float64}
    dict = Dict{NTuple{2,Float64},Float64}()
    for line in eachline(input_file)
        stripped_line = strip(line)
        isempty(stripped_line) && continue
        key_text, value_text = split_pi_line(stripped_line)
        dict[parse_pi_key(key_text)] = parse(Float64, strip(value_text))
    end
    return dict
end