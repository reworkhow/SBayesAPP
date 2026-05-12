module ConfigTypes 

export MarkerProbitTreeState, NonMPIConfig

const VALID_ANNOTATION_PRIOR_MODELS = (:group_dirichlet, :marker_probit_tree)

function normalize_annotation_prior_model(annotation_prior_model)
    model = annotation_prior_model isa Symbol ? annotation_prior_model : Symbol(annotation_prior_model)
    model in VALID_ANNOTATION_PRIOR_MODELS || error(
        "annotation_prior_model must be one of $(collect(VALID_ANNOTATION_PRIOR_MODELS)), got: $annotation_prior_model",
    )
    return model
end

struct NonMPIConfig
    data_path::String
    analysis_path::String
    nIter::Int
    seed::Int
    nrank::Int
    annot_file::String
    annot_dict::String
    out_freq::Int
    starting_value_dir::String
    gscale_value_dir::String
    st_path::String
    thin::Int
    n1::Int
    n2::Int
    n_con::Int
    annotation_prior_model::Symbol
    estimate_vare::Bool
    estimate_vara::Bool
    estimate_pi::Bool
    estimate_Gscale::Bool
    estGscale_iter::Int
    report_pleiotropic_qtl_effect_matrix::Bool
    output_mcmc_delta::Bool
    is_continue::Bool
end

mutable struct MarkerProbitTreeState
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
end

function NonMPIConfig(
    data_path::String,
    analysis_path::String,
    nIter::Int,
    seed::Int,
    nrank::Int,
    annot_file::String,
    annot_dict::String,
    out_freq::Int,
    starting_value_dir::String,
    gscale_value_dir::String,
    st_path::String,
    thin::Int,
    n1::Int,
    n2::Int,
    is_continue::Bool;
    n_con::Int=0,
    annotation_prior_model::Symbol=:group_dirichlet,
    estimate_vare::Bool=true,
    estimate_vara::Bool=true,
    estimate_pi::Bool=true,
    estimate_Gscale::Bool=true,
    estGscale_iter::Int=500,
    report_pleiotropic_qtl_effect_matrix::Bool=true,
    output_mcmc_delta::Bool=true,
)
    normalized_annotation_prior_model = normalize_annotation_prior_model(annotation_prior_model)
    return NonMPIConfig(
        data_path,
        analysis_path,
        nIter,
        seed,
        nrank,
        annot_file,
        annot_dict,
        out_freq,
        starting_value_dir,
        gscale_value_dir,
        st_path,
        thin,
        n1,
        n2,
        n_con,
        normalized_annotation_prior_model,
        estimate_vare,
        estimate_vara,
        estimate_pi,
        estimate_Gscale,
        estGscale_iter,
        report_pleiotropic_qtl_effect_matrix,
        output_mcmc_delta,
        is_continue,
    )
end

end