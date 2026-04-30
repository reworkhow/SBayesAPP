function initialize_r_blk_state(my_nblk, my_rank, Rprior; is_continue=false, starting_value_dir=nothing, fixed_r_dir=nothing)
    if !is_continue
        return [Rprior for _ in 1:my_nblk]
    end

    if !isnothing(fixed_r_dir)
        R_blk = [zeros(eltype(Rprior), size(Rprior)...) for _ in 1:my_nblk]
        for b in 1:my_nblk
            R_blk[b] = readdlm(fixed_r_dir * "estR.txt")
        end
        println("Loaded and Fixed R_blk from estR.txt in $fixed_r_dir for rank $my_rank")
        return R_blk
    end

    R_blk_folder = starting_value_dir * "last_sample_R_blk/"
    if isdir(R_blk_folder)
        R_blk = [zeros(eltype(Rprior), size(Rprior)...) for _ in 1:my_nblk]
        for b in 1:my_nblk
            R_blk[b] = readdlm(R_blk_folder * "R_blk_b$(b)_rank$(my_rank).txt")
        end
        println("Loaded R_blk from $R_blk_folder for rank $my_rank")
        return R_blk
    end

    println("Warning: Folder $R_blk_folder in rank $my_rank does not exist. R_blk will be initialized with Rprior.")
    return [Rprior for _ in 1:my_nblk]
end