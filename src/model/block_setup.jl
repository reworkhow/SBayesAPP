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

function build_all_marker_mask(anno_matrix_dict)
    return Dict(key => trues(size(anno_matrix_dict[key], 1), 1) for key in keys(anno_matrix_dict))
end

function build_global_annotation_design(blkID, blkSNPsIndex_dict, anno_matrix_dict; add_intercept::Bool=true)
    nmarker = sum(length(blkSNPsIndex_dict[blk]) for blk in blkID)
    nfeature = size(anno_matrix_dict[first(blkID)], 2)
    design_matrix = zeros(Float64, nmarker, nfeature + (add_intercept ? 1 : 0))

    if add_intercept
        design_matrix[:, 1] .= 1.0
    end

    for blk in blkID
        rows = blkSNPsIndex_dict[blk]
        block_annotations = Float64.(anno_matrix_dict[blk])
        if add_intercept
            design_matrix[rows, 2:end] .= block_annotations
        else
            design_matrix[rows, :] .= block_annotations
        end
    end

    return design_matrix
end

function prepare_block_state!(blkID, blkSNPsIndex_dict, transformed_x_dict, anno_matrix_dict, nCon, annotationType)
    reorder_block_snp_indices!(blkID, blkSNPsIndex_dict)
    xpx_dict, xArray_dict = build_block_designs(transformed_x_dict, anno_matrix_dict, nCon)
    anno_mask_dict = build_annotation_sampling_mask(anno_matrix_dict, nCon, annotationType)
    return xpx_dict, xArray_dict, anno_mask_dict
end

function prepare_marker_probit_tree_block_state!(blkID, blkSNPsIndex_dict, transformed_x_dict, anno_matrix_dict)
    reorder_block_snp_indices!(blkID, blkSNPsIndex_dict)
    xpx_dict, xArray_dict = build_block_designs(transformed_x_dict, anno_matrix_dict, 0)
    anno_mask_dict = build_all_marker_mask(anno_matrix_dict)
    annotation_design = build_global_annotation_design(blkID, blkSNPsIndex_dict, anno_matrix_dict; add_intercept=true)
    return xpx_dict, xArray_dict, anno_mask_dict, annotation_design
end