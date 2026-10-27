function build_nonmpi_sampler_settings(config::ConfigTypes.NonMPIConfig)
    uses_marker_tree = config.annotation_prior_model in (
        :marker_probit_tree,
        :marker_probit_tree_bayesr,
    )
    effective_n_con = config.n_con
    estimate_vara = config.estimate_vara
    estimate_pi = config.estimate_pi
    if uses_marker_tree && effective_n_con != 0
        @warn "n_con is ignored for marker-probit-tree priors; all annotation columns are treated as prior features."
        effective_n_con = 0
    end
    if uses_marker_tree && !estimate_pi
        @warn "estimate_pi=false is ignored for marker-probit-tree priors."
        estimate_pi = true
    end
    is_continue = config.is_continue
    if uses_marker_tree && is_continue
        @warn "is_continue=true is ignored for marker-probit-tree priors. A fresh run will be started instead."
        is_continue = false
    end
    return (
        estimate_vare=config.estimate_vare,
        estimate_vara=estimate_vara,
        estimate_pi=estimate_pi,
        estimate_Gscale=config.estimate_Gscale && estimate_vara,
        estGscale_iter=config.estGscale_iter,
        gscale_structure=config.gscale_structure,
        effective_n_con=effective_n_con,
        is_continue=is_continue,
    )
end
