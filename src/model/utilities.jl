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

function read_to_dict(input_file::String)::Dict{Vector{Float64},Float64}
    dict = Dict{Vector{Float64},Float64}()
    pi_index = 1
    for line in eachline(input_file)
        key_value = split(line, r"],", limit=2)

        if length(key_value) == 2
            value = parse(Float64, strip(key_value[2]))

            if pi_index == 1
                key = [0.0; 0.0]
            elseif pi_index == 2
                key = [1.0; 1.0]
            elseif pi_index == 3
                key = [1.0; 0.0]
            else
                key = [0.0; 1.0]
            end

            dict[key] = value
        end
        pi_index += 1
    end
    return dict
end

function read_to_dict_posterior_mean(input_file::String)::Dict{Vector{Float64},Float64}
    dict = Dict{Vector{Float64},Float64}()

    open(input_file, "r") do file
        for line in eachline(file)
            parts = split(line, '\t')
            key = eval(Meta.parse(parts[1]))
            value = parse(Float64, parts[2])
            dict[key] = value
        end
    end

    return dict
end