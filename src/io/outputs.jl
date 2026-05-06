using DelimitedFiles

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
)
    (; mean_pi, mean_pi2, meanB2, meanA2, meanBcor2, meanAcor2, meanA, meanAcor, meanB, meanBcor, meanG, meanG2, meanGcor, meanGcor2, meanGtotal, meanGtotal2, mcmcAtruecor_c, mcmcBcor_c, mcmcGcov_c, mcmcGcor_c, mcmcGcov_total, mcmcGcor_total, meanR, meanR2) = posterior_mean_state

    for cat in 1:nCategory
        writedlm(analysis_path * "estB" * string(cat) * ".txt", meanB[cat])
        meanBcor[cat] = mean(mcmcBcor_c[:, cat][.!isnan.(mcmcBcor_c[:, cat])])
        if report_pleiotropic_qtl_effect_matrix
            writedlm(analysis_path * "estA" * string(cat) * ".txt", meanA[cat])
            meanAcor[cat] = mean(mcmcAtruecor_c[:, cat][.!isnan.(mcmcAtruecor_c[:, cat])])
        end
        if estimate_vara
            writedlm(analysis_path * "estB_std" * string(cat) * ".txt", sqrt.((meanB2[cat] .- (meanB[cat] .^ 2))))
            meanBcor2[cat] = mean(mcmcBcor_c[:, cat][.!isnan.(mcmcBcor_c[:, cat])] .^ 2)
            if report_pleiotropic_qtl_effect_matrix
                writedlm(analysis_path * "estA_std" * string(cat) * ".txt", sqrt.((meanA2[cat] .- (meanA[cat] .^ 2))))
                meanAcor2[cat] = mean(mcmcAtruecor_c[:, cat][.!isnan.(mcmcAtruecor_c[:, cat])] .^ 2)
            end
        end
    end

    writedlm(analysis_path * "mcmcGcov_c.txt", mcmcGcov_c)
    writedlm(analysis_path * "mcmcGcov_total.txt", mcmcGcov_total)
    writedlm(analysis_path * "mcmcGcor_c.txt", mcmcGcor_c)
    writedlm(analysis_path * "mcmcGcor_total.txt", mcmcGcor_total)

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

    for cat in 1:nCategory
        writedlm(analysis_path * "estG" * string(cat) * ".txt", meanG[cat])
        writedlm(analysis_path * "estG_std" * string(cat) * ".txt", sqrt.(abs.(meanG2[cat] .- (meanG[cat] .^ 2))))
        meanGcor[cat] = mean(mcmcGcor_c[:, cat][.!isnan.(mcmcGcor_c[:, cat])])
        meanGcor2[cat] = mean(mcmcGcor_c[:, cat][.!isnan.(mcmcGcor_c[:, cat])] .^ 2)
    end
    writedlm(analysis_path * "estGcor.txt", meanGcor)
    writedlm(analysis_path * "estGcor_std.txt", sqrt.(meanGcor2 .- (meanGcor .^ 2)))

    writedlm(analysis_path * "estGtotal.txt", meanGtotal, ',')
    writedlm(analysis_path * "estGtotal_std.txt", sqrt.((meanGtotal2 .- (meanGtotal .^ 2))), ',')
    meanGcor_total = mean(mcmcGcor_total[.!isnan.(mcmcGcor_total)])
    meanGcor_total2 = mean(mcmcGcor_total[.!isnan.(mcmcGcor_total)] .^ 2)
    writedlm(analysis_path * "estGcor_total.txt", meanGcor_total)
    writedlm(analysis_path * "estGcor_total_std.txt", sqrt(meanGcor_total2 - (meanGcor_total^2)))

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

    return nothing
end