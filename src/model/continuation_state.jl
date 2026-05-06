function load_effect_state(effect_starting_path, delta_starting_path, my_rank, nTraits)
    betaArray = [
        vec(readdlm(effect_starting_path * "last_mcmc_betaArray$(trait).rank$my_rank.txt"))
        for trait in 1:nTraits
    ]
    deltaArray = [
        vec(readdlm(delta_starting_path * "last_sample_delta$(trait)_rank$my_rank.txt"))
        for trait in 1:nTraits
    ]
    alphaArray = [deltaArray[trait] .* betaArray[trait] for trait in 1:nTraits]
    return betaArray, alphaArray, deltaArray
end

function apply_alpha_correction!(
    transformed_y_dict,
    xArray_dict,
    blkID,
    blkSNPsIndex_dict,
    alphaArray,
    my_nsnp,
    nCategory,
    nTraits,
)
    for blk in blkID
        base_x = xArray_dict[blk][end]
        nMarkerb = size(base_x, 2)
        for trait in 1:nTraits
            alphaArray_total = zeros(nMarkerb)
            for category in 1:nCategory
                alphaArray_total += alphaArray[trait][(category - 1) * my_nsnp .+ blkSNPsIndex_dict[blk]]
            end
            transformed_y_dict[blk][trait] = transformed_y_dict[blk][trait] - base_x * alphaArray_total
        end
    end
    return nothing
end

function initialize_effect_state!(
    effect_starting_path,
    delta_starting_path,
    my_rank,
    my_nsnp,
    nCategory,
    nTraits,
    blkID,
    blkSNPsIndex_dict,
    xArray_dict,
    transformed_y_dict,
)
    betaArray, alphaArray, deltaArray = load_effect_state(
        effect_starting_path,
        delta_starting_path,
        my_rank,
        nTraits,
    )
    apply_alpha_correction!(
        transformed_y_dict,
        xArray_dict,
        blkID,
        blkSNPsIndex_dict,
        alphaArray,
        my_nsnp,
        nCategory,
        nTraits,
    )
    return betaArray, alphaArray, deltaArray
end