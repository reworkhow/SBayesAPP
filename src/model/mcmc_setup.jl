function zero_pi_moments(Pi)
    mean_pi = deepcopy(Pi)
    mean_pi2 = deepcopy(Pi)
    for category in eachindex(mean_pi)
        for key in keys(mean_pi[category])
            mean_pi[category][key] = 0.0
            mean_pi2[category][key] = 0.0
        end
    end
    return mean_pi, mean_pi2
end

function initialize_rank0_mcmc_state(
    Pi,
    nIter,
    burnin,
    thin,
    nTraits,
    nCategory;
    estimate_pi=true,
    estimate_vara=true,
    estimate_vare=true,
    report_pleiotropic_qtl_effect_matrix=true,
    save_category_correlation_outputs=true,
)
    nsample4mean = Int(floor((nIter - burnin) / thin))
    mean_pi = nothing
    mean_pi2 = nothing
    if estimate_pi
        mean_pi, mean_pi2 = zero_pi_moments(Pi)
    end

    meanB2 = estimate_vara ? [zeros(nTraits, nTraits) for _ in 1:nCategory] : nothing
    meanA2 = estimate_vara && report_pleiotropic_qtl_effect_matrix ? [zeros(nTraits, nTraits) for _ in 1:nCategory] : nothing
    meanBcor2 = estimate_vara && save_category_correlation_outputs ? zeros(nCategory) : nothing
    meanAcor2 = estimate_vara && report_pleiotropic_qtl_effect_matrix && save_category_correlation_outputs ? zeros(nCategory) : nothing

    meanA = report_pleiotropic_qtl_effect_matrix ? [zeros(nTraits, nTraits) for _ in 1:nCategory] : nothing
    meanAcor = report_pleiotropic_qtl_effect_matrix && save_category_correlation_outputs ? zeros(nCategory) : nothing
    meanB = [zeros(nTraits, nTraits) for _ in 1:nCategory]
    meanBcor = save_category_correlation_outputs ? zeros(nCategory) : nothing
    meanG = [zeros(nTraits, nTraits) for _ in 1:nCategory]
    meanG2 = [zeros(nTraits, nTraits) for _ in 1:nCategory]
    meanGcor = save_category_correlation_outputs ? zeros(nCategory) : nothing
    meanGcor2 = save_category_correlation_outputs ? zeros(nCategory) : nothing
    meanGcor_count = save_category_correlation_outputs ? zeros(Int, nCategory) : nothing
    meanSSE = [zeros(nTraits, nTraits) for _ in 1:nCategory]
    meanGtotal = zeros(nTraits, nTraits)
    meanGtotal2 = zeros(nTraits, nTraits)
    meanGcor_total = 0.0
    meanGcor_total2 = 0.0
    meanGcor_total_count = 0
    mcmcAtruecor_c = report_pleiotropic_qtl_effect_matrix && save_category_correlation_outputs ? zeros(nsample4mean, nCategory) : nothing
    mcmcBcor_c = save_category_correlation_outputs ? zeros(nsample4mean, nCategory) : nothing
    meanR = estimate_vare ? zeros(nTraits, nTraits) : nothing
    meanR2 = estimate_vare ? zeros(nTraits, nTraits) : nothing

    return (
        nsample4mean=nsample4mean,
        mean_pi=mean_pi,
        mean_pi2=mean_pi2,
        meanB2=meanB2,
        meanA2=meanA2,
        meanBcor2=meanBcor2,
        meanAcor2=meanAcor2,
        meanA=meanA,
        meanAcor=meanAcor,
        meanB=meanB,
        meanBcor=meanBcor,
        meanG=meanG,
        meanG2=meanG2,
        meanGcor=meanGcor,
        meanGcor2=meanGcor2,
        meanGcor_count=meanGcor_count,
        meanSSE=meanSSE,
        meanGtotal=meanGtotal,
        meanGtotal2=meanGtotal2,
        meanGcor_total=meanGcor_total,
        meanGcor_total2=meanGcor_total2,
        meanGcor_total_count=meanGcor_total_count,
        mcmcAtruecor_c=mcmcAtruecor_c,
        mcmcBcor_c=mcmcBcor_c,
        meanR=meanR,
        meanR2=meanR2,
    )
end

function prepare_mcmc_output_files(analysis_path; report_pleiotropic_qtl_effect_matrix=true)
    file_names = Dict(
        "pi" => analysis_path * "MCMC_samples_pi.txt",
        "beta_effects_variance" => analysis_path * "MCMC_samples_beta_effects_variance.txt",
        "genetic_effects_variance" => analysis_path * "MCMC_samples_genetic_effects_variance.txt",
        "total_genetic_effects_variance" => analysis_path * "MCMC_samples_total_genetic_effects_variance.txt",
    )
    if report_pleiotropic_qtl_effect_matrix
        file_names["marker_effects_variance"] = analysis_path * "MCMC_samples_marker_effects_variance.txt"
    end
    for path in values(file_names)
        if isfile(path)
            println("File $path already exists! It will be overwritten.")
        else
            println("Creating file: $path")
        end
        open(path, "w") do io
        end
    end
    return file_names
end