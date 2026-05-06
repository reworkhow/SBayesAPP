using LinearAlgebra: cholesky, diag

function flatten(vec_of_vecs)
    return vcat(map(x -> x[:], vec_of_vecs)...)
end

function unflatten(seq_vec, subvec_length)
    return [seq_vec[i:i+subvec_length-1] for i in 1:subvec_length:length(seq_vec)]
end

function flatten_matrices(vec_of_mats)
    return vcat(map(x -> x[:], vec_of_mats)...)
end

function unflatten_matrices(seq_vec, nrows, ncols)
    num_matrices = length(seq_vec) ÷ (nrows * ncols)
    return [reshape(seq_vec[(i-1)*nrows*ncols+1:i*nrows*ncols], nrows, ncols) for i in 1:num_matrices]
end

function compute_correlation(cov_matrix)
    std_devs = sqrt.(diag(cov_matrix))
    return cov_matrix[1, 2] / (std_devs[1] * std_devs[2])
end

function is_positive_definite(matrix)
    try
        cholesky(matrix)
        return true
    catch
        return false
    end
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

read_to_dict_posterior_mean(input_file::String)::Dict{NTuple{2,Float64},Float64} = read_to_dict(input_file)