function build_nonmpi_sampler_settings(config::ConfigTypes.NonMPIConfig)
    effective_n_con = config.n_con
    estimate_vara = config.estimate_vara
    estimate_pi = config.estimate_pi
    if config.annotation_prior_model == :marker_probit_tree && effective_n_con != 0
        @warn "n_con is ignored when annotation_prior_model=:marker_probit_tree; all annotation columns are treated as prior features."
        effective_n_con = 0
    end
    if config.annotation_prior_model == :marker_probit_tree && !estimate_pi
        @warn "estimate_pi=false is ignored when annotation_prior_model=:marker_probit_tree."
        estimate_pi = true
    end
    is_continue = config.is_continue
    if config.annotation_prior_model == :marker_probit_tree && is_continue
        @warn "is_continue=true is ignored when annotation_prior_model=:marker_probit_tree. A fresh run will be started instead."
        is_continue = false
    end
    return (
        estimate_vare=config.estimate_vare,
        estimate_vara=estimate_vara,
        estimate_pi=estimate_pi,
        estimate_Gscale=config.estimate_Gscale && estimate_vara,
        estGscale_iter=config.estGscale_iter,
        effective_n_con=effective_n_con,
        is_continue=is_continue,
    )
end