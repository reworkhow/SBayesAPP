using DelimitedFiles: writedlm

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
    mcmc_Delta,
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
        last_delta = mcmc_Delta[trait][:, end]
        open(analysis_path * "last_sample_delta/last_sample_delta$(trait)_rank$(my_rank).txt", "w") do io
            writedlm(io, last_delta)
        end
    end
end