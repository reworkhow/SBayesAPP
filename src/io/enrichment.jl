using Statistics: mean, std

const VALID_ENRICHMENT_PARAMETERS = (:coh2, :h21, :h22)
const PI_STATE_KEYS = ("pi00", "pi11", "pi10", "pi01")

function normalize_enrichment_parameter(parameter)
    normalized = parameter isa Symbol ? parameter : Symbol(lowercase(String(parameter)))
    normalized in VALID_ENRICHMENT_PARAMETERS || error(
        "parameter must be one of $(collect(VALID_ENRICHMENT_PARAMETERS)), got: $parameter",
    )
    return normalized
end


function extract_genetic_matrix_parameter(matrix_values::AbstractMatrix{<:Real}, parameter::Symbol)
    parameter === :coh2 && return Float64(matrix_values[1, 2])
    parameter === :h21 && return Float64(matrix_values[1, 1])
    parameter === :h22 && return Float64(matrix_values[2, 2])
    error("Unsupported enrichment parameter: $parameter")
end

function read_matrix_series(path::AbstractString)
    isfile(path) || error("Missing MCMC file: $path")
    values = readdlm(path, ',')
    matrix_values = ndims(values) == 1 ? reshape(values, :, 1) : values
    size(matrix_values, 2) == 2 || error("Expected 2 columns in $path, got $(size(matrix_values, 2))")
    size(matrix_values, 1) % 2 == 0 || error("Expected an even number of rows in $path")
    return Float64.(matrix_values)
end

function read_category_parameter_samples(path::AbstractString, n_annotation::Int, parameter::Symbol)
    matrix_values = read_matrix_series(path)
    n_matrices = div(size(matrix_values, 1), 2)
    n_matrices % n_annotation == 0 || error(
        "File $path contains $n_matrices matrices, which is not divisible by n_annotation=$n_annotation",
    )

    n_samples = div(n_matrices, n_annotation)
    values = Matrix{Float64}(undef, n_samples, n_annotation)

    for sample_index in 1:n_samples
        for annotation_index in 1:n_annotation
            matrix_index = (sample_index - 1) * n_annotation + annotation_index
            row_start = (matrix_index - 1) * 2 + 1
            values[sample_index, annotation_index] = extract_genetic_matrix_parameter(
                @view(matrix_values[row_start:(row_start + 1), :]),
                parameter,
            )
        end
    end

    return values
end

function read_total_parameter_samples(path::AbstractString, parameter::Symbol)
    matrix_values = read_matrix_series(path)
    n_samples = div(size(matrix_values, 1), 2)
    values = Vector{Float64}(undef, n_samples)

    for sample_index in 1:n_samples
        row_start = (sample_index - 1) * 2 + 1
        values[sample_index] = extract_genetic_matrix_parameter(
            @view(matrix_values[row_start:(row_start + 1), :]),
            parameter,
        )
    end

    return values
end

function extract_genetic_matrix_correlation(matrix_values::AbstractMatrix{<:Real})
    varg1 = Float64(matrix_values[1, 1])
    varg2 = Float64(matrix_values[2, 2])
    gcov = Float64(matrix_values[1, 2])
    return (varg1 > 0.0 && varg2 > 0.0) ? gcov / sqrt(varg1 * varg2) : NaN
end

function read_category_correlation_samples(path::AbstractString, n_annotation::Int)
    matrix_values = read_matrix_series(path)
    n_matrices = div(size(matrix_values, 1), 2)
    n_matrices % n_annotation == 0 || error(
        "File $path contains $n_matrices matrices, which is not divisible by n_annotation=$n_annotation",
    )

    n_samples = div(n_matrices, n_annotation)
    values = Matrix{Float64}(undef, n_samples, n_annotation)

    for sample_index in 1:n_samples
        for annotation_index in 1:n_annotation
            matrix_index = (sample_index - 1) * n_annotation + annotation_index
            row_start = (matrix_index - 1) * 2 + 1
            values[sample_index, annotation_index] = extract_genetic_matrix_correlation(
                @view(matrix_values[row_start:(row_start + 1), :]),
            )
        end
    end

    return values
end

function read_category_pi_samples(path::AbstractString, n_annotation::Int)
    isfile(path) || error("Missing MCMC file: $path")
    lines = collect(eachline(path))
    n_state = length(PI_STATE_KEYS)
    n_record_lines = length(lines)
    n_record_lines % (n_annotation * n_state) == 0 || error(
        "File $path contains $n_record_lines Pi lines, which is not divisible by n_annotation * $n_state",
    )

    n_samples = div(n_record_lines, n_annotation * n_state)
    values = Array{Float64}(undef, n_samples, n_annotation, n_state)
    expected_keys = pi_key_order()

    for sample_index in 1:n_samples
        for annotation_index in 1:n_annotation
            for state_index in 1:n_state
                line_index = ((sample_index - 1) * n_annotation + (annotation_index - 1)) * n_state + state_index
                key_text, value_text = split_pi_line(strip(lines[line_index]))
                parsed_key = parse_pi_key(key_text)
                parsed_key == expected_keys[state_index] || error(
                    "Unexpected Pi state order in $path at line $line_index: expected $(format_pi_key(expected_keys[state_index])), got $key_text",
                )
                values[sample_index, annotation_index, state_index] = parse(Float64, strip(value_text))
            end
        end
    end

    return values
end

function resolve_marker_effects_variance_file(chain_path::AbstractString)
    marker_file = joinpath(chain_path, "MCMC_samples_marker_effects_variance.txt")
    isfile(marker_file) && return marker_file

    beta_file = joinpath(chain_path, "MCMC_samples_beta_effects_variance.txt")
    isfile(beta_file) && return beta_file

    error(
        "Missing marker/beta effect variance MCMC file in $chain_path. Expected MCMC_samples_marker_effects_variance.txt or MCMC_samples_beta_effects_variance.txt",
    )
end

function load_annotation_enrichment_inputs(data_path::AbstractString, annot_file::AbstractString)
    annot_path = resolve_input_path(data_path, annot_file)
    annot = read_input_table(annot_path)
    "SNP" in names(annot) || error("Annotation file must contain a SNP column: $annot_path")

    annotation_columns = names(annot)[2:end]
    isempty(annotation_columns) && error("Annotation file must contain at least one annotation column: $annot_path")

    annotation_values = Float64.(Matrix(annot[:, annotation_columns]))
    row_weight_totals = vec(sum(annotation_values, dims=2))
    any(iszero, row_weight_totals) && error(
        "Annotation file contains SNP rows with zero total annotation weight: $annot_path",
    )

    per_snp_normalization_weights = 1.0 ./ row_weight_totals
    weighted_annotation_counts = vec(transpose(annotation_values) * per_snp_normalization_weights)

    return (
        annotation_path=annot_path,
        annotation_names=String.(annotation_columns),
        annotation_table=annot,
        weighted_annotation_counts=weighted_annotation_counts,
        n_total_effects=nrow(annot),
    )
end

function classify_enrichment(enrichment_pp::Real, depletion_pp::Real, significance_threshold::Real)
    enrichment_pp > significance_threshold && return "Enrichment"
    depletion_pp > significance_threshold && return "Depletion"
    return "Non-significant"
end

function build_single_chain_enrichment_tables(
    category_values,
    total_values,
    annotation_names,
    annotation_table,
    weighted_annotation_counts,
    significance_threshold,
)
    conditional_values = vec(mean(category_values, dims=1))
    conditional_posterior_gcov = DataFrame(annotationName=annotation_names, conditional_gcov=conditional_values)
    sort!(conditional_posterior_gcov, :conditional_gcov, rev=true)

    n_total_effects = nrow(annotation_table)
    relative_values = category_values ./ reshape(total_values, :, 1)
    posterior_enrichment_summary = DataFrame(
        Annotation=annotation_names,
        Mean=zeros(length(annotation_names)),
        SD=zeros(length(annotation_names)),
        EnrPP=zeros(length(annotation_names)),
        DeplPP=zeros(length(annotation_names)),
        Significance=fill("", length(annotation_names)),
    )

    for annotation_index in eachindex(annotation_names)
        annotation_enrichment_samples =
            relative_values[:, annotation_index] .* (n_total_effects / weighted_annotation_counts[annotation_index])
        posterior_enrichment_summary.Mean[annotation_index] = mean(annotation_enrichment_samples)
        posterior_enrichment_summary.SD[annotation_index] = std(annotation_enrichment_samples)
        posterior_enrichment_summary.EnrPP[annotation_index] = count(abs.(annotation_enrichment_samples) .> 1.0) / length(annotation_enrichment_samples)
        posterior_enrichment_summary.DeplPP[annotation_index] = count(abs.(annotation_enrichment_samples) .< 1.0) / length(annotation_enrichment_samples)
        posterior_enrichment_summary.Significance[annotation_index] = classify_enrichment(
            posterior_enrichment_summary.EnrPP[annotation_index],
            posterior_enrichment_summary.DeplPP[annotation_index],
            significance_threshold,
        )
    end

    sort!(posterior_enrichment_summary, :Mean, rev=true)
    return (
        conditional_posterior_gcov=conditional_posterior_gcov,
        posterior_enrichment_summary=posterior_enrichment_summary,
    )
end

function parse_chain_paths(chain_paths::Union{AbstractString,AbstractVector{<:AbstractString}})
    paths = chain_paths isa AbstractString ? split(chain_paths, ',') : collect(chain_paths)
    normalized = [String(strip(path)) for path in paths if !isempty(strip(path))]
    isempty(normalized) && error("At least one chain path is required")
    for path in normalized
        isdir(path) || error("Chain path not found: $path")
    end
    return normalized
end

function discover_chain_paths(group_path::AbstractString; chain_prefix::AbstractString="seed_")
    isdir(group_path) || error("Group path not found: $group_path")
    entries = filter(name -> startswith(name, chain_prefix), readdir(group_path))
    isempty(entries) && error("No chain directories starting with '$chain_prefix' found in $group_path")
    sort!(entries, by=name -> begin
        suffix = replace(name, chain_prefix => ""; count=1)
        parsed = tryparse(Int, suffix)
        parsed === nothing ? typemax(Int) : parsed
    end)
    return [joinpath(group_path, name) for name in entries]
end

posterior_sample_indices(n_samples::Int) = collect(1:n_samples)

function assert_same_sample_count(sample_counts, file_label::AbstractString)
    first_count = first(sample_counts)
    all(==(first_count), sample_counts) || error(
        "Sample count mismatch across chains for $file_label: $(join(sample_counts, ", "))",
    )
    return first_count
end

function group_average_matrices(matrices::Vector{Matrix{Float64}}, indices::AbstractVector{Int})
    sample_counts = [size(matrix, 1) for matrix in matrices]
    assert_same_sample_count(sample_counts, "matrix samples")
    feature_counts = [size(matrix, 2) for matrix in matrices]
    all(==(first(feature_counts)), feature_counts) || error(
        "Feature count mismatch across chains: $(join(feature_counts, ", "))",
    )

    averaged = zeros(Float64, length(indices), first(feature_counts))
    for matrix in matrices
        averaged .+= @view(matrix[indices, :])
    end
    averaged ./= length(matrices)
    return averaged
end

function group_average_matrices_skip_nonfinite(matrices::Vector{Matrix{Float64}}, indices::AbstractVector{Int})
    sample_counts = [size(matrix, 1) for matrix in matrices]
    assert_same_sample_count(sample_counts, "matrix samples")
    feature_counts = [size(matrix, 2) for matrix in matrices]
    all(==(first(feature_counts)), feature_counts) || error(
        "Feature count mismatch across chains: $(join(feature_counts, ", "))",
    )

    averaged = fill(NaN, length(indices), first(feature_counts))
    for iteration_index in eachindex(indices)
        sample_index = indices[iteration_index]
        for feature_index in 1:first(feature_counts)
            finite_values = [matrix[sample_index, feature_index] for matrix in matrices if isfinite(matrix[sample_index, feature_index])]
            if !isempty(finite_values)
                averaged[iteration_index, feature_index] = mean(finite_values)
            end
        end
    end
    return averaged
end

function group_average_vectors(vectors::Vector{Vector{Float64}}, indices::AbstractVector{Int})
    sample_counts = [length(vector) for vector in vectors]
    assert_same_sample_count(sample_counts, "vector samples")

    averaged = zeros(Float64, length(indices))
    for vector in vectors
        averaged .+= @view(vector[indices])
    end
    averaged ./= length(vectors)
    return averaged
end

function group_average_pi_arrays(arrays::Vector{Array{Float64,3}}, indices::AbstractVector{Int})
    sample_counts = [size(array, 1) for array in arrays]
    assert_same_sample_count(sample_counts, "Pi samples")
    annotation_counts = [size(array, 2) for array in arrays]
    all(==(first(annotation_counts)), annotation_counts) || error(
        "Pi annotation count mismatch across chains: $(join(annotation_counts, ", "))",
    )

    averaged = zeros(Float64, length(indices), first(annotation_counts), length(PI_STATE_KEYS))
    for array in arrays
        averaged .+= @view(array[indices, :, :])
    end
    averaged ./= length(arrays)
    return averaged
end

function write_iteration_table(path::AbstractString, posterior_samples, values::AbstractMatrix, column_names)
    mkpath(dirname(path))
    table = DataFrame(PosteriorSample=posterior_samples)
    for (index, column_name) in enumerate(column_names)
        table[!, Symbol(column_name)] = values[:, index]
    end
    CSV.write(path, table)
    return path
end

function build_parameter_summary(values::AbstractMatrix, annotation_names; mean_name::AbstractString="Mean", sd_name::AbstractString="SD", skip_nonfinite::Bool=false)
    summary = DataFrame(Annotation=annotation_names)
    means = Vector{Float64}(undef, length(annotation_names))
    sds = Vector{Float64}(undef, length(annotation_names))
    for annotation_index in eachindex(annotation_names)
        samples = values[:, annotation_index]
        if skip_nonfinite
            samples = samples[isfinite.(samples)]
        end
        if isempty(samples)
            means[annotation_index] = NaN
            sds[annotation_index] = NaN
        else
            means[annotation_index] = mean(samples)
            sds[annotation_index] = std(samples)
        end
    end
    summary[!, Symbol(mean_name)] = means
    summary[!, Symbol(sd_name)] = sds
    sort!(summary, Symbol(mean_name), rev=true)
    return summary
end

function build_pi_summary(pi_values::Array{Float64,3}, annotation_names)
    summary = DataFrame(Annotation=String[])
    for state_key in PI_STATE_KEYS
        summary[!, Symbol("$(state_key)_mean")] = Float64[]
        summary[!, Symbol("$(state_key)_sd")] = Float64[]
    end
    for column_name in (
        "trait1_polygenicity_mean",
        "trait1_polygenicity_sd",
        "trait2_polygenicity_mean",
        "trait2_polygenicity_sd",
        "pi11_standardized_mean",
        "pi11_standardized_sd",
        "pi10_standardized_mean",
        "pi10_standardized_sd",
        "pi01_standardized_mean",
        "pi01_standardized_sd",
    )
        summary[!, Symbol(column_name)] = Float64[]
    end

    for annotation_index in eachindex(annotation_names)
        pi00 = pi_values[:, annotation_index, 1]
        pi11 = pi_values[:, annotation_index, 2]
        pi10 = pi_values[:, annotation_index, 3]
        pi01 = pi_values[:, annotation_index, 4]
        trait1 = pi11 .+ pi10
        trait2 = pi11 .+ pi01
        any_trait = 1 .- pi00
        pi11_standardized = pi11 ./ any_trait
        pi10_standardized = pi10 ./ any_trait
        pi01_standardized = pi01 ./ any_trait

        push!(summary, (
            annotation_names[annotation_index],
            mean(pi00), std(pi00),
            mean(pi11), std(pi11),
            mean(pi10), std(pi10),
            mean(pi01), std(pi01),
            mean(trait1), std(trait1),
            mean(trait2), std(trait2),
            mean(pi11_standardized), std(pi11_standardized),
            mean(pi10_standardized), std(pi10_standardized),
            mean(pi01_standardized), std(pi01_standardized),
        ))
    end

    return summary
end

default_multi_chain_output_dir(group_path::AbstractString) = joinpath(group_path, "multi_chain_summary")

function calculate_multi_chain_enrichment(
    data_path::AbstractString,
    annot_file::AbstractString;
    group_path::Union{Nothing,AbstractString}=nothing,
    chain_paths::Union{Nothing,AbstractString,AbstractVector{<:AbstractString}}=nothing,
    output_path::Union{Nothing,AbstractString}=nothing,
    chain_prefix::AbstractString="seed_",
    group_label::AbstractString="Group1",
    enrichment_parameter::Union{Symbol,AbstractString}=:coh2,
    abs_total::Bool=false,
    significance_threshold::Real=0.90,
    write_iteration_averages::Bool=true,
)
    0.0 <= significance_threshold <= 1.0 || error(
        "PP significance_threshold must be between 0 and 1, got: $significance_threshold",
    )
    normalized_parameter = normalize_enrichment_parameter(enrichment_parameter)
    resolved_chain_paths = if chain_paths !== nothing
        parse_chain_paths(chain_paths)
    elseif group_path !== nothing
        discover_chain_paths(String(group_path); chain_prefix=chain_prefix)
    else
        error("Either group_path or chain_paths must be provided")
    end
    resolved_group_path = group_path === nothing ? dirname(first(resolved_chain_paths)) : String(group_path)
    resolved_output_path = output_path === nothing ? default_multi_chain_output_dir(resolved_group_path) : String(output_path)
    isabspath(resolved_output_path) || (resolved_output_path = joinpath(resolved_group_path, resolved_output_path))
    mkpath(resolved_output_path)

    annotation_inputs = load_annotation_enrichment_inputs(data_path, annot_file)
    annotation_names = annotation_inputs.annotation_names
    n_annotation = length(annotation_names)

    category_matrices = [
        read_category_parameter_samples(
            joinpath(chain_path, "MCMC_samples_genetic_effects_variance.txt"),
            n_annotation,
            normalized_parameter,
        ) for chain_path in resolved_chain_paths
    ]
    total_vectors = [
        read_total_parameter_samples(
            joinpath(chain_path, "MCMC_samples_total_genetic_effects_variance.txt"),
            normalized_parameter,
        ) for chain_path in resolved_chain_paths
    ]

    n_samples = assert_same_sample_count([size(matrix, 1) for matrix in category_matrices], "genetic-variance samples")
    assert_same_sample_count([length(vector) for vector in total_vectors], "total genetic-variance samples")
    sample_indices = posterior_sample_indices(n_samples)

    # Average chains at each posterior sample before summarizing across samples.
    category_group_average = group_average_matrices(category_matrices, sample_indices)
    total_group_average = group_average_vectors(total_vectors, sample_indices)
    abs_total && (total_group_average = abs.(total_group_average))

    enrichment_tables = build_single_chain_enrichment_tables(
        category_group_average,
        total_group_average,
        annotation_names,
        annotation_inputs.annotation_table,
        annotation_inputs.weighted_annotation_counts,
        significance_threshold,
    )

    snpcor_matrices = [
        read_category_correlation_samples(
            resolve_marker_effects_variance_file(chain_path),
            n_annotation,
        ) for chain_path in resolved_chain_paths
    ]
    # SNP correlations are undefined when a variance is nonpositive; match the old R SNPcor workflow by skipping those values.
    snpcor_group_average = group_average_matrices_skip_nonfinite(snpcor_matrices, sample_indices)
    snpcor_summary = build_parameter_summary(snpcor_group_average, annotation_names; mean_name="snpcor_mean", sd_name="snpcor_sd", skip_nonfinite=true)

    pi_arrays = [
        read_category_pi_samples(joinpath(chain_path, "MCMC_samples_pi.txt"), n_annotation)
        for chain_path in resolved_chain_paths
    ]
    pi_group_average = group_average_pi_arrays(pi_arrays, sample_indices)
    pi_summary = build_pi_summary(pi_group_average, annotation_names)

    prefix = "$(group_label)_$(length(resolved_chain_paths))chains"
    enrichment_summary_file = joinpath(resolved_output_path, "$(prefix)_$(normalized_parameter)_enrichment_summary.csv")
    conditional_file = joinpath(resolved_output_path, "$(prefix)_conditional_$(normalized_parameter).csv")
    snpcor_file = joinpath(resolved_output_path, "$(prefix)_snpcor_summary.csv")
    pi_file = joinpath(resolved_output_path, "$(prefix)_pi_summary.csv")

    CSV.write(enrichment_summary_file, enrichment_tables.posterior_enrichment_summary)
    CSV.write(conditional_file, enrichment_tables.conditional_posterior_gcov)
    CSV.write(snpcor_file, snpcor_summary)
    CSV.write(pi_file, pi_summary)

    average_files = String[]
    if write_iteration_averages
        push!(
            average_files,
            write_iteration_table(
                joinpath(resolved_output_path, "$(prefix)_avg_$(normalized_parameter)_MCMC.csv"),
                sample_indices,
                category_group_average,
                annotation_names,
            ),
        )
        push!(
            average_files,
            write_iteration_table(
                joinpath(resolved_output_path, "$(prefix)_avg_snpcor_MCMC.csv"),
                sample_indices,
                snpcor_group_average,
                annotation_names,
            ),
        )
        pi_average_table = DataFrame(PosteriorSample=sample_indices)
        for (annotation_index, annotation_name) in enumerate(annotation_names)
            for (state_index, state_key) in enumerate(PI_STATE_KEYS)
                pi_average_table[!, Symbol("$(annotation_name)_$(state_key)")] = pi_group_average[:, annotation_index, state_index]
            end
        end
        pi_average_file = joinpath(resolved_output_path, "$(prefix)_avg_pi_MCMC.csv")
        CSV.write(pi_average_file, pi_average_table)
        push!(average_files, pi_average_file)
    end

    return (
        group_path=resolved_group_path,
        chain_paths=resolved_chain_paths,
        output_path=resolved_output_path,
        group_label=String(group_label),
        n_chains=length(resolved_chain_paths),
        nsamples=n_samples,
        enrichment_parameter=normalized_parameter,
        enrichment_summary=enrichment_tables.posterior_enrichment_summary,
        conditional_summary=enrichment_tables.conditional_posterior_gcov,
        snpcor_summary=snpcor_summary,
        pi_summary=pi_summary,
        enrichment_summary_file=enrichment_summary_file,
        conditional_summary_file=conditional_file,
        snpcor_summary_file=snpcor_file,
        pi_summary_file=pi_file,
        average_files=average_files,
    )
end

function default_single_chain_summary_filename(parameter::Symbol)
    return "single_chain_$(parameter)_summary.csv"
end

function default_single_chain_conditional_filename(parameter::Symbol)
    return "conditional_$(parameter).csv"
end

function calculate_single_chain_enrichment(
    data_path::AbstractString,
    analysis_path::AbstractString,
    annot_file::AbstractString;
    parameter::Union{Symbol,AbstractString}=:coh2,
    abs_total::Bool=false,
    significance_threshold::Real=0.90,
    summary_output_file::Union{Nothing,AbstractString}=nothing,
    conditional_output_file::Union{Nothing,AbstractString}=nothing,
)
    normalized_parameter = normalize_enrichment_parameter(parameter)
    analysis_dir = String(analysis_path)
    isdir(analysis_dir) || error("Analysis path not found: $analysis_dir")
    0.0 <= significance_threshold <= 1.0 || error(
        "PP significance_threshold must be between 0 and 1, got: $significance_threshold",
    )

    annotation_inputs = load_annotation_enrichment_inputs(data_path, annot_file)
    n_annotation = length(annotation_inputs.annotation_names)

    category_file = joinpath(analysis_dir, "MCMC_samples_genetic_effects_variance.txt")
    total_file = joinpath(analysis_dir, "MCMC_samples_total_genetic_effects_variance.txt")
    category_values = read_category_parameter_samples(category_file, n_annotation, normalized_parameter)
    total_values = read_total_parameter_samples(total_file, normalized_parameter)
    size(category_values, 1) == length(total_values) || error(
        "Sample count mismatch between $category_file and $total_file",
    )

    abs_total && (total_values = abs.(total_values))

    tables = build_single_chain_enrichment_tables(
        category_values,
        total_values,
        annotation_inputs.annotation_names,
        annotation_inputs.annotation_table,
        annotation_inputs.weighted_annotation_counts,
        significance_threshold,
    )

    summary_path = if summary_output_file === nothing
        joinpath(analysis_dir, default_single_chain_summary_filename(normalized_parameter))
    else
        summary_candidate = String(summary_output_file)
        isabspath(summary_candidate) ? summary_candidate : joinpath(analysis_dir, summary_candidate)
    end
    conditional_path = if conditional_output_file === nothing
        joinpath(analysis_dir, default_single_chain_conditional_filename(normalized_parameter))
    else
        conditional_candidate = String(conditional_output_file)
        isabspath(conditional_candidate) ? conditional_candidate : joinpath(analysis_dir, conditional_candidate)
    end

    mkpath(dirname(summary_path))
    mkpath(dirname(conditional_path))
    CSV.write(summary_path, tables.posterior_enrichment_summary)
    CSV.write(conditional_path, tables.conditional_posterior_gcov)

    return (
        analysis_path=analysis_dir,
        annotation_path=annotation_inputs.annotation_path,
        parameter=normalized_parameter,
        nsamples=size(category_values, 1),
        significance_threshold=Float64(significance_threshold),
        posterior_enrichment_summary=tables.posterior_enrichment_summary,
        conditional_posterior_gcov=tables.conditional_posterior_gcov,
        posterior_enrichment_summary_file=summary_path,
        conditional_posterior_gcov_file=conditional_path,
    )
end