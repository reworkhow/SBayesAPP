using Statistics: mean, std

const VALID_SINGLE_CHAIN_ENRICHMENT_PARAMETERS = (:coh2, :h21, :h22)

function normalize_single_chain_enrichment_parameter(parameter)
    normalized = parameter isa Symbol ? parameter : Symbol(lowercase(String(parameter)))
    normalized in VALID_SINGLE_CHAIN_ENRICHMENT_PARAMETERS || error(
        "parameter must be one of $(collect(VALID_SINGLE_CHAIN_ENRICHMENT_PARAMETERS)), got: $parameter",
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
    normalized_parameter = normalize_single_chain_enrichment_parameter(parameter)
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