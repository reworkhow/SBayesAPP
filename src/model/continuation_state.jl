function apply_alpha_correction!(
    blocks,
    alphaArray,
    my_nsnp,
    nCategory,
    nTraits,
)
    for block in blocks
        base_x = block.x_arrays[end]
        nMarkerb = size(base_x, 2)
        for trait in 1:nTraits
            alphaArray_total = zeros(nMarkerb)
            for category in 1:nCategory
                alphaArray_total += alphaArray[trait][(category - 1) * my_nsnp .+ block.snp_indices]
            end
            block.transformed_y[trait] = block.transformed_y[trait] - base_x * alphaArray_total
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
    blocks,
)
    betaArray, alphaArray, deltaArray = load_effect_state(
        effect_starting_path,
        delta_starting_path,
        my_rank,
        nTraits,
    )
    apply_alpha_correction!(
        blocks,
        alphaArray,
        my_nsnp,
        nCategory,
        nTraits,
    )
    return betaArray, alphaArray, deltaArray
end