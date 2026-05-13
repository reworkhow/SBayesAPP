function build_block_designs(block::NonMPIBlockData, nCon)
    Xb = block.transformed_x
    nMarkerb = size(Xb, 2)
    xpx = Vector{Vector{Float64}}(undef, nCon + 1)
    x_arrays = Vector{Matrix{Float64}}(undef, nCon + 1)

    for c in 1:nCon
        annot_weights = block.anno_matrix[:, c]
        xpx[c] = [annot_weights[i]^2 * dot(Xb[:, i], Xb[:, i]) for i in 1:nMarkerb]
        x_arrays[c] = Xb * diagm(annot_weights)
    end

    xpx[end] = [dot(Xb[:, i], Xb[:, i]) for i in 1:nMarkerb]
    x_arrays[end] = Xb
    return xpx, x_arrays
end

function build_annotation_sampling_mask(block::NonMPIBlockData, nCon, annotationType)
    annotation_mask = block.anno_matrix .!= 0.0
    if nCon > 0
        continuous_columns = findall(annotationType .== "continue")
        annotation_mask[:, continuous_columns] .= true
    end
    return annotation_mask
end

build_all_marker_mask(block::NonMPIBlockData) = trues(size(block.anno_matrix, 1), 1)

function build_global_annotation_design(blocks; add_intercept::Bool=true)
    nmarker = sum(length(block.snp_indices) for block in blocks)
    nfeature = size(first(blocks).anno_matrix, 2)
    design_matrix = zeros(Float64, nmarker, nfeature + (add_intercept ? 1 : 0))

    if add_intercept
        design_matrix[:, 1] .= 1.0
    end

    for block in blocks
        rows = block.snp_indices
        block_annotations = Float64.(block.anno_matrix)
        if add_intercept
            design_matrix[rows, 2:end] .= block_annotations
        else
            design_matrix[rows, :] .= block_annotations
        end
    end

    return design_matrix
end

function prepare_block_state!(blocks, nCon, annotationType)
    for block in blocks
        block.xpx, block.x_arrays = build_block_designs(block, nCon)
        block.annotation_mask = build_annotation_sampling_mask(block, nCon, annotationType)
    end
    return nothing
end

function prepare_marker_probit_tree_block_state!(blocks)
    for block in blocks
        block.xpx, block.x_arrays = build_block_designs(block, 0)
        block.annotation_mask = build_all_marker_mask(block)
    end
    annotation_design = build_global_annotation_design(blocks; add_intercept=true)
    return annotation_design
end