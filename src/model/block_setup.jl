function reorder_block_snp_indices!(blkID, blkSNPsIndex_dict)
    sort!(blkID)
    nblk = length(blkID)
    for i in 2:nblk
        previous_block = blkID[i - 1]
        cumulative_snp_count = length(blkSNPsIndex_dict[previous_block])
        for key in blkID[i:nblk]
            blkSNPsIndex_dict[key] = blkSNPsIndex_dict[key] .+ cumulative_snp_count
        end
    end
    return nothing
end

function build_block_designs(transformed_x_dict, anno_matrix_dict, nCon)
    xpx_dict = Dict{Int,Vector{Vector{Float64}}}()
    xArray_dict = Dict{Int,Vector{Matrix{Float64}}}()

    for blk in keys(transformed_x_dict)
        Xb = transformed_x_dict[blk]
        nMarkerb = size(Xb, 2)
        xpx_dict[blk] = Vector{Vector{Float64}}(undef, nCon + 1)
        xArray_dict[blk] = Vector{Matrix{Float64}}(undef, nCon + 1)

        for c in 1:nCon
            annot_weights = anno_matrix_dict[blk][:, c]
            xpx_dict[blk][c] = [annot_weights[i]^2 * dot(Xb[:, i], Xb[:, i]) for i in 1:nMarkerb]
            xArray_dict[blk][c] = Xb * diagm(annot_weights)
        end

        xpx_dict[blk][end] = [dot(Xb[:, i], Xb[:, i]) for i in 1:nMarkerb]
        xArray_dict[blk][end] = Xb
    end

    return xpx_dict, xArray_dict
end

function build_annotation_sampling_mask(anno_matrix_dict, nCon, annotationType)
    anno_mask_dict = Dict(key => (anno_matrix_dict[key] .!= 0.0) for key in keys(anno_matrix_dict))
    if nCon > 0
        continuous_columns = findall(annotationType .== "continue")
        for blk in keys(anno_mask_dict)
            anno_mask_dict[blk][:, continuous_columns] .= true
        end
    end
    return anno_mask_dict
end

function prepare_block_state!(blkID, blkSNPsIndex_dict, transformed_x_dict, anno_matrix_dict, nCon, annotationType)
    reorder_block_snp_indices!(blkID, blkSNPsIndex_dict)
    xpx_dict, xArray_dict = build_block_designs(transformed_x_dict, anno_matrix_dict, nCon)
    anno_mask_dict = build_annotation_sampling_mask(anno_matrix_dict, nCon, annotationType)
    return xpx_dict, xArray_dict, anno_mask_dict
end