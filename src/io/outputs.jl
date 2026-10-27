using DelimitedFiles

function write_nonmpi_restart_state(
    analysis_path,
    my_rank,
    my_nblk,
    nCategory,
    nTraits,
    betaArray,
    R_blk,
    A_vec,
    Pi,
    deltaArray,
)
    for trait in 1:nTraits
        writedlm(analysis_path * "last_mcmc_betaArray$trait.rank$my_rank.txt", betaArray[trait])
    end

    mkpath(analysis_path * "last_sample_R_blk/")
    for block in 1:my_nblk
        writedlm(analysis_path * "last_sample_R_blk/R_blk_b$(block)_rank$(my_rank).txt", R_blk[block])
    end

    mkpath(analysis_path * "beta_effect_var_matrices_last_sample/")
    for category in 1:nCategory
        writedlm(analysis_path * "beta_effect_var_matrices_last_sample/beta_effect_matrix_$(category).txt", A_vec[category], ',')
    end

    mkpath(analysis_path * "pi_last_sample/")
    for category in 1:nCategory
        write_pi_dict(analysis_path * "pi_last_sample/pi_$(category).txt", Pi[category])
    end

    mkpath(analysis_path * "last_sample_delta/")
    for trait in 1:nTraits
        writedlm(analysis_path * "last_sample_delta/last_sample_delta$(trait)_rank$(my_rank).txt", deltaArray[trait])
    end
end

function write_annotation_probit_coefficients(output_file::AbstractString, annotation_names, step_labels, values)
    open(output_file, "w") do io
        println(io, "Annotation\tStep\tEstimate")
        for (annotation_index, annotation_name) in enumerate(annotation_names)
            for (step_index, step_label) in enumerate(step_labels)
                println(io, annotation_name, '\t', step_label, '\t', values[annotation_index, step_index])
            end
        end
    end
    return nothing
end

function write_bayesr_joint_prior_summary(output_file::AbstractString, means, standard_deviations)
    open(output_file, "w") do io
        println(io, "Trait1Class\tTrait2Class\tEstimate\tSD")
        for state in marker_probit_tree_bayesr_state_order()
            index = marker_probit_tree_bayesr_state_index(state)
            println(io, state[1], '\t', state[2], '\t', means[index], '\t', standard_deviations[index])
        end
    end
    return nothing
end

function sanitize_gene_output_component(value::AbstractString)
    component = replace(strip(String(value)), r"[^A-Za-z0-9._-]+" => "_")
    component = strip(component, ['_', '.', '-'])
    return isempty(component) ? "category" : component
end

function write_gene_covariance_outputs(analysis_path::AbstractString, state::GeneCovarianceState)
    output_dir = normalize_dir_path(analysis_path)
    category_component = sanitize_gene_output_component(state.category_name)
    weight_component = String(state.mapping_weight)
    output_stem = "gene_genetic_covariance_$(category_component)_$(weight_component)"
    main_path = output_dir * output_stem * ".txt"
    qc_path = output_dir * output_stem * "_qc.txt"

    open(main_path, "w") do io
        println(
            io,
            "GENE\tN_CATEGORY_X_SNPS\tEFFECTIVE_SNP_COUNT\tPOSTERIOR_MEAN_COV\tPOSTERIOR_SD_COV",
        )
        for gene_index in eachindex(state.gene_names)
            posterior_mean = state.posterior_mean_cov[gene_index]
            posterior_sd = sqrt(abs(
                state.posterior_second_moment_cov[gene_index] - posterior_mean^2,
            ))
            println(
                io,
                state.gene_names[gene_index],
                '\t',
                state.snp_counts[gene_index],
                '\t',
                state.effective_snp_counts[gene_index],
                '\t',
                posterior_mean,
                '\t',
                posterior_sd,
            )
        end
    end

    qc = state.qc
    qc_entries = [
        "RESOLVED_CATEGORY_NAME" => state.category_name,
        "RESOLVED_CATEGORY_INDEX" => state.category_index,
        "GENE_MAPPING_WEIGHT" => String(state.mapping_weight),
        "VALID_UNIQUE_SNP_GENE_PAIRS" => qc.valid_unique_snp_gene_pairs,
        "VALID_CATEGORY_X_MAPPED_SNPS" => qc.valid_category_mapped_snps,
        "CATEGORY_X_SNPS_WITHOUT_VALID_GENE_MAPPING" => qc.category_snps_without_valid_gene_mapping,
        "SNPS_MAPPED_TO_ONE_RETAINED_GENE" => qc.snps_mapped_to_one_retained_gene,
        "SNPS_MAPPED_TO_MULTIPLE_RETAINED_GENES" => qc.snps_mapped_to_multiple_retained_genes,
        "MAXIMUM_RETAINED_GENES_PER_SNP" => qc.maximum_retained_genes_per_snp,
        "RETAINED_GENES" => qc.retained_genes,
        "RETAINED_MCMC_SAMPLES" => state.retained_sample_count,
    ]
    open(qc_path, "w") do io
        println(io, "METRIC\tVALUE")
        for (metric, value) in qc_entries
            println(io, metric, '\t', value)
        end
    end
    return (main=main_path, qc=qc_path)
end

function save_nonmpi_posterior_mean!(
    analysis_path,
    my_rank,
    nCategory,
    meanAlpha,
    posterior_mean_state;
    estimate_vara,
    estimate_vare,
    estimate_pi,
    report_pleiotropic_qtl_effect_matrix,
    annotation_prior_model=:group_dirichlet,
    marker_probit_tree_state=nothing,
    annotation_names=String[],
)
    (; mean_pi, mean_pi2, meanB2, meanA2, meanBcor2, meanAcor2, meanA, meanAcor, meanB, meanBcor, meanG, meanG2, meanGcor, meanGcor2, meanGcor_count, meanGtotal, meanGtotal2, meanGcor_total, meanGcor_total2, meanGcor_total_count, mcmcAtruecor_c, mcmcBcor_c, meanR, meanR2) = posterior_mean_state
    save_category_correlation_outputs = annotation_prior_model == :group_dirichlet

    for cat in 1:nCategory
        writedlm(analysis_path * "estB" * string(cat) * ".txt", meanB[cat])
        if save_category_correlation_outputs
            meanBcor[cat] = mean(mcmcBcor_c[:, cat][.!isnan.(mcmcBcor_c[:, cat])])
        end
        if report_pleiotropic_qtl_effect_matrix
            writedlm(analysis_path * "estA" * string(cat) * ".txt", meanA[cat])
            if save_category_correlation_outputs
                meanAcor[cat] = mean(mcmcAtruecor_c[:, cat][.!isnan.(mcmcAtruecor_c[:, cat])])
            end
        end
        if estimate_vara
            writedlm(analysis_path * "estB_std" * string(cat) * ".txt", sqrt.((meanB2[cat] .- (meanB[cat] .^ 2))))
            if save_category_correlation_outputs
                meanBcor2[cat] = mean(mcmcBcor_c[:, cat][.!isnan.(mcmcBcor_c[:, cat])] .^ 2)
            end
            if report_pleiotropic_qtl_effect_matrix
                writedlm(analysis_path * "estA_std" * string(cat) * ".txt", sqrt.((meanA2[cat] .- (meanA[cat] .^ 2))))
                if save_category_correlation_outputs
                    meanAcor2[cat] = mean(mcmcAtruecor_c[:, cat][.!isnan.(mcmcAtruecor_c[:, cat])] .^ 2)
                end
            end
        end
    end

    if save_category_correlation_outputs
        writedlm(analysis_path * "estBcor.txt", meanBcor)
        if report_pleiotropic_qtl_effect_matrix
            writedlm(analysis_path * "mcmcAtruecor_c.txt", mcmcAtruecor_c)
            writedlm(analysis_path * "estAcor.txt", meanAcor)
        end

        if estimate_vara
            writedlm(analysis_path * "estBcor_std.txt", sqrt.(abs.(meanBcor2 .- (meanBcor .^ 2))))
            if report_pleiotropic_qtl_effect_matrix
                writedlm(analysis_path * "estAcor_std.txt", sqrt.(abs.(meanAcor2 .- (meanAcor .^ 2))))
            end
        end
    end

    for cat in 1:nCategory
        writedlm(analysis_path * "estG" * string(cat) * ".txt", meanG[cat])
        writedlm(analysis_path * "estG_std" * string(cat) * ".txt", sqrt.(abs.(meanG2[cat] .- (meanG[cat] .^ 2))))
    end
    if save_category_correlation_outputs
        for cat in 1:nCategory
            if meanGcor_count[cat] > 0
                meanGcor[cat] /= meanGcor_count[cat]
                meanGcor2[cat] /= meanGcor_count[cat]
            else
                meanGcor[cat] = NaN
                meanGcor2[cat] = NaN
            end
        end
        writedlm(analysis_path * "estGcor.txt", meanGcor)
        writedlm(analysis_path * "estGcor_std.txt", sqrt.(abs.(meanGcor2 .- (meanGcor .^ 2))))
    end

    writedlm(analysis_path * "estGtotal.txt", meanGtotal, ',')
    writedlm(analysis_path * "estGtotal_std.txt", sqrt.((meanGtotal2 .- (meanGtotal .^ 2))), ',')
    if meanGcor_total_count > 0
        meanGcor_total /= meanGcor_total_count
        meanGcor_total2 /= meanGcor_total_count
    else
        meanGcor_total = NaN
        meanGcor_total2 = NaN
    end
    writedlm(analysis_path * "estGcor_total.txt", meanGcor_total)
    writedlm(analysis_path * "estGcor_total_std.txt", sqrt(abs(meanGcor_total2 - (meanGcor_total^2))))

    if estimate_vare
        writedlm(analysis_path * "estR.txt", meanR)
        writedlm(analysis_path * "estR_std.txt", sqrt.((meanR2 .- (meanR .^ 2))))
    end

    if estimate_pi
        for cat in 1:nCategory
            write_pi_dict(analysis_path * "estPi" * string(cat) * ".txt", mean_pi[cat])
            std_pi = deepcopy(mean_pi[cat])
            for key in pi_key_order()
                std_pi[key] = sqrt(mean_pi2[cat][key] - mean_pi[cat][key]^2)
            end
            write_pi_dict(analysis_path * "estPi_std" * string(cat) * ".txt", std_pi)
        end
    end

    for trait in eachindex(meanAlpha)
        writedlm(analysis_path * "meanAlpha" * string(trait) * ".rank$my_rank.txt", meanAlpha[trait])
    end

    if is_marker_probit_tree_model(annotation_prior_model) && marker_probit_tree_state !== nothing
        coefficient_names = isempty(annotation_names) ?
            ["Intercept"; ["Annotation_$index" for index in 1:(size(marker_probit_tree_state.mean_coefficients, 1) - 1)]] :
            ["Intercept"; annotation_names]
        coefficient_sd = sqrt.(abs.(marker_probit_tree_state.mean_coefficients2 .- (marker_probit_tree_state.mean_coefficients .^ 2)))
        write_annotation_probit_coefficients(
            analysis_path * "annotation_probit_coefficients.txt",
            coefficient_names,
            marker_probit_tree_step_labels(marker_probit_tree_state),
            marker_probit_tree_state.mean_coefficients,
        )
        write_annotation_probit_coefficients(
            analysis_path * "annotation_probit_coefficients_std.txt",
            coefficient_names,
            marker_probit_tree_step_labels(marker_probit_tree_state),
            coefficient_sd,
        )
        if is_marker_probit_tree_bayesr(annotation_prior_model)
            prior_sd = sqrt.(abs.(
                marker_probit_tree_state.mean_snp_pi2 .-
                marker_probit_tree_state.mean_snp_pi .^ 2
            ))
            write_bayesr_joint_prior_summary(
                analysis_path * "annotation_bayesr_predicted_joint_state_probabilities.txt",
                marker_probit_tree_state.mean_snp_pi,
                prior_sd,
            )
        end
    end

    return nothing
end
