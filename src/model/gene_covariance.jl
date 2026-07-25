struct GeneBlockSlice
    gene_index::Int
    marker_indices::Vector{Int}
    effect_indices::Vector{Int}
    weights::Vector{Float64}
end

struct GeneCovarianceQC
    valid_unique_snp_gene_pairs::Int
    valid_category_mapped_snps::Int
    category_snps_without_valid_gene_mapping::Int
    snps_mapped_to_one_retained_gene::Int
    snps_mapped_to_multiple_retained_genes::Int
    maximum_retained_genes_per_snp::Int
    retained_genes::Int
end

mutable struct GeneCovarianceState
    gene_names::Vector{String}
    category_name::String
    category_index::Int
    mapping_weight::Symbol
    per_block_slices::Vector{Vector{GeneBlockSlice}}
    snp_counts::Vector{Int}
    effective_snp_counts::Vector{Float64}
    posterior_mean_cov::Vector{Float64}
    posterior_second_moment_cov::Vector{Float64}
    retained_sample_count::Int
    qc::GeneCovarianceQC
end

function resolve_gene_category(annotation_metadata, category_name::AbstractString)
    requested_name = String(category_name)
    matches = findall(==(requested_name), annotation_metadata.annotationName)
    length(matches) == 1 || error(
        isempty(matches) ?
        "Gene category '$requested_name' is not an exact annotation name." :
        "Gene category '$requested_name' is ambiguous because it occurs more than once.",
    )

    category_index = only(matches)
    annotation_metadata.annotationType[category_index] == "category" || error(
        "Gene category '$requested_name' must be categorical; continuous annotations are not supported.",
    )
    return category_index
end

function validate_binary_gene_category(blocks, category_index::Int, category_name::AbstractString)
    for block in blocks
        values = view(block.anno_matrix, :, category_index)
        for (marker_index, value) in enumerate(values)
            isfinite(value) || error(
                "Gene category '$category_name' contains a non-finite value in block $(block.block_id), marker $marker_index.",
            )
            (value == 0.0 || value == 1.0) || error(
                "Gene category '$category_name' must be binary (0/1); found $value in block $(block.block_id), marker $marker_index.",
            )
        end
    end
    return nothing
end

function read_snp_gene_map(mapping_path::AbstractString)
    delim, ignore_repeated = detect_input_delimiter(mapping_path)
    mapping = CSV.read(mapping_path, DataFrame; delim=delim, ignorerepeated=ignore_repeated)
    all(column -> column in names(mapping), ("SNP", "GENE")) || error(
        "SNP-to-gene mapping must contain columns named SNP and GENE.",
    )

    pairs = Tuple{String,String}[]
    for (row_index, row) in enumerate(eachrow(mapping))
        raw_snp = row.SNP
        raw_gene = row.GENE
        (ismissing(raw_snp) || ismissing(raw_gene)) && continue
        snp = strip(string(raw_snp))
        gene = strip(string(raw_gene))
        (isempty(snp) || isempty(gene)) && continue
        (occursin(r"[\t\r\n]", snp) || occursin(r"[\t\r\n]", gene)) && error(
            "SNP-to-gene mapping row $row_index contains a tab or newline inside SNP or GENE.",
        )
        push!(pairs, (snp, gene))
    end

    return unique(pairs)
end

function build_gene_snp_locations(blocks, snp_ids_dict)
    locations = Dict{String,Tuple{Int,Int}}()
    for (block_position, block) in enumerate(blocks)
        block_id = Int(block.block_id)
        snp_ids = snp_ids_dict[block_id]
        for (marker_index, snp_id) in enumerate(snp_ids)
            snp = String(snp_id)
            haskey(locations, snp) && error(
                "Analysis SNP ID '$snp' occurs more than once across LD blocks; gene mapping would be ambiguous.",
            )
            locations[snp] = (block_position, marker_index)
        end
    end
    return locations
end

# creates a dictionary that tells, for a given block and category, which compact effect index corresponds to each SNP/marker index.
function category_effect_index_lookup(block, category_index::Int)
    marker_indices = block.category_marker_indices[category_index]
    effect_indices = block.category_effect_indices[category_index]
    length(marker_indices) == length(effect_indices) || error(
        "Internal category marker/effect index mismatch for block $(block.block_id), category $category_index.",
    )
    lookup = Dict{Int,Int}()
    for position in eachindex(marker_indices)
        lookup[marker_indices[position]] = effect_indices[position]
    end
    return lookup
end

function initialize_gene_covariance_state(
    blocks,
    snp_ids_dict,
    annotation_metadata,
    mapping_path::AbstractString,
    category_name::AbstractString,
    mapping_weight::Symbol,
)
    category_index = resolve_gene_category(annotation_metadata, category_name)
    validate_binary_gene_category(blocks, category_index, category_name)

    mapping_pairs = read_snp_gene_map(mapping_path)
    locations = build_gene_snp_locations(blocks, snp_ids_dict)
    category_membership = Dict{String,Bool}()
    category_snps = Set{String}()
    for (snp, (block_position, marker_index)) in locations
        is_member = blocks[block_position].anno_matrix[marker_index, category_index] != 0.0
        category_membership[snp] = is_member
        is_member && push!(category_snps, snp)
    end

    # keeps only SNP-gene pairs whose SNP is present in the analysis data and belongs to the chosen category.
    valid_pairs = [
        pair for pair in mapping_pairs
        if haskey(locations, pair[1]) && category_membership[pair[1]]
    ]

    gene_names = sort!(collect(Set(last(pair) for pair in valid_pairs)))
    gene_indices = Dict(gene => index for (index, gene) in enumerate(gene_names))
    # each SNP maps to the set of genes it belongs to.
    genes_by_snp = Dict{String,Set{String}}()
    for (snp, gene) in valid_pairs
        push!(get!(genes_by_snp, snp, Set{String}()), gene)
    end
    # counts how many genes each SNP maps to.
    gene_multiplicity = Dict(snp => length(genes) for (snp, genes) in genes_by_snp)

    snps_by_gene = [Set{String}() for _ in gene_names]
    # keyed by: (block_position, gene_index)
    # values: vector of entries (marker_index, effect_index, weight)
    entries_by_slice = Dict{Tuple{Int,Int},Vector{Tuple{Int,Int,Float64}}}()
    effect_lookups = [
        category_effect_index_lookup(block, category_index)
        for block in blocks
    ]
    for (snp, gene) in valid_pairs
        block_position, marker_index = locations[snp]
        gene_index = gene_indices[gene]
        push!(snps_by_gene[gene_index], snp)
        effect_index = get(effect_lookups[block_position], marker_index, 0)
        effect_index != 0 || error(
            "Internal error: category-X SNP '$snp' has no compact category effect index.",
        )
        # may need derivations later
        weight = mapping_weight == :full ? 1.0 : inv(sqrt(gene_multiplicity[snp]))
        push!(
            get!(entries_by_slice, (block_position, gene_index), Tuple{Int,Int,Float64}[]),
            (marker_index, effect_index, weight),
        )
    end

    # converts the grouped SNP-gene entries into per-block GeneBlockSlice objects.
    per_block_slices = [GeneBlockSlice[] for _ in blocks]
    for ((block_position, gene_index), entries) in entries_by_slice
        sort!(entries; by=first)
        push!(
            per_block_slices[block_position],
            GeneBlockSlice(
                gene_index, # which gene this slice belongs to
                first.(entries), # all marker indices for that gene in this block.
                getindex.(entries, 2), # all effect indices.
                getindex.(entries, 3), # all SNP-gene mapping weights.
            ),
        )
    end 
    # sorts the gene slices inside each block by gene_index
    for slices in per_block_slices
        sort!(slices; by=slice -> slice.gene_index)
    end

    # computes per-gene SNP counts and QC summary statistics for the SNP-to-gene mapping.
    snp_counts = length.(snps_by_gene)
    # computes a weighted SNP count per gene.
    effective_snp_counts = [
        sum(
            mapping_weight == :full ? 1.0 : inv(gene_multiplicity[snp])
            for snp in gene_snps
        )
        for gene_snps in snps_by_gene
    ]
    valid_category_mapped_snps = Set(keys(genes_by_snp))
    multiplicities = collect(values(gene_multiplicity))
    qc = GeneCovarianceQC(
        length(valid_pairs), # unique pairs retained after SNP and category filtering
        length(valid_category_mapped_snps), # analysis category-X SNPs mapped to retained genes
        length(setdiff(category_snps, valid_category_mapped_snps)), # category-X SNPs without a valid gene
        count(==(1), multiplicities), # SNPs mapped to exactly one retained gene
        count(>=(2), multiplicities), # SNPs mapped to multiple retained genes
        isempty(multiplicities) ? 0 : maximum(multiplicities), # largest retained-gene count for one SNP
        length(gene_names), # genes containing at least one valid category-X SNP
    )

    return GeneCovarianceState(
        gene_names,
        String(category_name),
        category_index,
        mapping_weight,
        per_block_slices,
        snp_counts,
        effective_snp_counts,
        zeros(length(gene_names)),
        zeros(length(gene_names)),
        0,
        qc,
    )
end

function accumulate_gene_covariance_sample!(state::GeneCovarianceState, blocks, alphaArray)
    sample_covariance = zeros(length(state.gene_names))
    category_index = state.category_index

    for (block_position, slices) in enumerate(state.per_block_slices)
        isempty(slices) && continue
        transformed_x = blocks[block_position].transformed_x
        n_eigen = size(transformed_x, 1)
        for slice in slices
            predictor_trait1 = zeros(n_eigen)
            predictor_trait2 = zeros(n_eigen)
            for position in eachindex(slice.marker_indices)
                marker_index = slice.marker_indices[position]
                effect_index = slice.effect_indices[position]
                effect_weight = slice.weights[position]
                alpha_trait1 = effect_weight * alphaArray[1][category_index][effect_index]
                alpha_trait2 = effect_weight * alphaArray[2][category_index][effect_index]
                alpha_trait1 != 0.0 && LinearAlgebra.axpy!(
                    alpha_trait1,
                    view(transformed_x, :, marker_index),
                    predictor_trait1,
                )
                alpha_trait2 != 0.0 && LinearAlgebra.axpy!(
                    alpha_trait2,
                    view(transformed_x, :, marker_index),
                    predictor_trait2,
                )
            end
            sample_covariance[slice.gene_index] += dot(predictor_trait1, predictor_trait2)
        end
    end

    state.retained_sample_count += 1
    moment_weight = inv(state.retained_sample_count)
    for gene_index in eachindex(sample_covariance)
        covariance = sample_covariance[gene_index]
        state.posterior_mean_cov[gene_index] +=
            (covariance - state.posterior_mean_cov[gene_index]) * moment_weight
        state.posterior_second_moment_cov[gene_index] +=
            (covariance^2 - state.posterior_second_moment_cov[gene_index]) * moment_weight
    end
    return sample_covariance
end
