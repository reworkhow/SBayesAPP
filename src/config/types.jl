module ConfigTypes 

export MarkerProbitTreeState, NonMPIConfig

const VALID_ANNOTATION_PRIOR_MODELS = (
    :group_dirichlet,
    :marker_probit_tree,
    :marker_probit_tree_bayesr,
)
const VALID_GENE_MAPPING_WEIGHTS = (:full, :fractional)
const VALID_GSCALE_STRUCTURES = (:correlated, :independent)

function normalize_annotation_prior_model(annotation_prior_model)
    model = annotation_prior_model isa Symbol ? annotation_prior_model : Symbol(annotation_prior_model)
    model in VALID_ANNOTATION_PRIOR_MODELS || error(
        "annotation_prior_model must be one of $(collect(VALID_ANNOTATION_PRIOR_MODELS)), got: $annotation_prior_model",
    )
    return model
end

function normalize_gene_mapping_weight(gene_mapping_weight)
    weight = gene_mapping_weight isa Symbol ? gene_mapping_weight : Symbol(lowercase(String(gene_mapping_weight)))
    weight in VALID_GENE_MAPPING_WEIGHTS || error(
        "gene_mapping_weight must be one of $(collect(VALID_GENE_MAPPING_WEIGHTS)), got: $gene_mapping_weight",
    )
    return weight
end

function normalize_gscale_structure(gscale_structure)
    structure = gscale_structure isa Symbol ? gscale_structure : Symbol(lowercase(String(gscale_structure)))
    structure in VALID_GSCALE_STRUCTURES || error(
        "gscale_structure must be correlated or independent, got: $gscale_structure",
    )
    return structure
end

function validate_nonmpi_schedule(nIter::Int, burnin::Int, thin::Int)
    nIter > 0 || error("nIter must be positive, got: $nIter")
    thin > 0 || error("thin must be positive, got: $thin")
    burnin >= 0 || error("burnin must be non-negative, got: $burnin")
    burnin < nIter || error("burnin must be smaller than nIter, got burnin=$burnin and nIter=$nIter")
    floor(Int, (nIter - burnin) / thin) >= 1 || error(
        "thin and burnin leave no post-burnin samples, got nIter=$nIter, burnin=$burnin, thin=$thin",
    )
    return nothing
end

struct NonMPIConfig
    data_path::String
    analysis_path::String
    nIter::Int
    seed::Int
    annot_file::String
    annot_dict::String
    starting_value_dir::String
    gscale_value_dir::String
    st_path::String
    thin::Int
    burnin::Int
    n1::Int
    n2::Int
    n_con::Int
    annotation_prior_model::Symbol
    estimate_vare::Bool
    estimate_vara::Bool
    estimate_pi::Bool
    estimate_Gscale::Bool
    estGscale_iter::Int
    gscale_structure::Symbol
    report_pleiotropic_qtl_effect_matrix::Bool
    output_mcmc_delta::Bool
    is_continue::Bool
    snp_gene_map::Union{Nothing,String}
    gene_category::Union{Nothing,String}
    gene_mapping_weight::Symbol
end

mutable struct MarkerProbitTreeState
    prior_model::Symbol
    design_matrix::Matrix{Float64}
    coefficients::Matrix{Float64}
    mean_coefficients::Matrix{Float64}
    mean_coefficients2::Matrix{Float64}
    variance::Vector{Float64}
    liability::Matrix{Float64}
    mu::Matrix{Float64}
    lower_bound::Matrix{Float64}
    upper_bound::Matrix{Float64}
    snp_pi::Matrix{Float64}
    mean_snp_pi::Vector{Float64}
    mean_snp_pi2::Vector{Float64}
end

function NonMPIConfig(
    data_path::String,
    analysis_path::String,
    nIter::Int,
    seed::Int,
    annot_file::String,
    annot_dict::String,
    starting_value_dir::String,
    gscale_value_dir::String,
    st_path::String,
    thin::Int,
    n1::Int,
    n2::Int,
    is_continue::Bool;
    burnin::Union{Nothing, Int}=nothing,
    n_con::Int=0,
    annotation_prior_model::Symbol=:group_dirichlet,
    estimate_vare::Bool=true,
    estimate_vara::Bool=true,
    estimate_pi::Bool=true,
    estimate_Gscale::Bool=true,
    estGscale_iter::Int=500,
    gscale_structure=:correlated,
    report_pleiotropic_qtl_effect_matrix::Bool=true,
    output_mcmc_delta::Bool=true,
    snp_gene_map::Union{Nothing,AbstractString}=nothing,
    gene_category::Union{Nothing,AbstractString}=nothing,
    gene_mapping_weight=:full,
)
    normalized_annotation_prior_model = normalize_annotation_prior_model(annotation_prior_model)
    normalized_gene_mapping_weight = normalize_gene_mapping_weight(gene_mapping_weight)
    normalized_gscale_structure = normalize_gscale_structure(gscale_structure)
    normalized_snp_gene_map = isnothing(snp_gene_map) ? nothing : String(snp_gene_map)
    normalized_gene_category = isnothing(gene_category) ? nothing : String(gene_category)
    has_snp_gene_map = !isnothing(normalized_snp_gene_map)
    has_gene_category = !isnothing(normalized_gene_category)
    has_snp_gene_map == has_gene_category || error(
        "--snp-gene-map and --gene-category must be provided together.",
    )
    if has_snp_gene_map
        isempty(strip(normalized_snp_gene_map)) && error("--snp-gene-map must not be empty.")
        isempty(strip(normalized_gene_category)) && error("--gene-category must not be empty.")
        normalized_annotation_prior_model == :group_dirichlet || error(
            "Gene-level covariance requires annotation_prior_model=:group_dirichlet.",
        )
        is_continue && error(
            "Gene-level covariance does not support is_continue=true in this MVP; restart integration is not implemented.",
        )
    elseif normalized_gene_mapping_weight != :full
        error("--gene-mapping-weight requires --snp-gene-map and --gene-category.")
    end
    resolved_burnin = isnothing(burnin) ? (is_continue ? 0 : floor(Int, nIter * 0.4)) : burnin
    validate_nonmpi_schedule(nIter, resolved_burnin, thin)
    return NonMPIConfig(
        data_path,
        analysis_path,
        nIter,
        seed,
        annot_file,
        annot_dict,
        starting_value_dir,
        gscale_value_dir,
        st_path,
        thin,
        resolved_burnin,
        n1,
        n2,
        n_con,
        normalized_annotation_prior_model,
        estimate_vare,
        estimate_vara,
        estimate_pi,
        estimate_Gscale,
        estGscale_iter,
        normalized_gscale_structure,
        report_pleiotropic_qtl_effect_matrix,
        output_mcmc_delta,
        is_continue,
        normalized_snp_gene_map,
        normalized_gene_category,
        normalized_gene_mapping_weight,
    )
end

end
