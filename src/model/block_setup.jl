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

function build_annotation_category_entries(annotation_mask::AbstractMatrix{Bool}, annotation_values::AbstractMatrix)
    nmarker, ncategory = size(annotation_mask)
    size(annotation_values) == size(annotation_mask) || error("annotation values must have the same shape as annotation_mask")
    category_indices = Vector{Vector{Int}}(undef, nmarker)
    category_values = Vector{Vector{Float64}}(undef, nmarker)
    for marker in 1:nmarker
        marker_categories = Int[]
        marker_values = Float64[]
        for category in 1:ncategory
            if annotation_mask[marker, category]
                push!(marker_categories, category)
                push!(marker_values, Float64(annotation_values[marker, category]))
            end
        end
        category_indices[marker] = marker_categories
        category_values[marker] = marker_values
    end
    return category_indices, category_values
end

function assign_annotation_effect_indices!(blocks, nCategory::Int)
    category_counts = zeros(Int, nCategory)

    for block in blocks
        marker_category_indices = block.annotation_category_indices::Vector{Vector{Int}}
        marker_effect_indices = Vector{Vector{Int}}(undef, length(marker_category_indices))
        block_category_marker_indices = [Int[] for _ in 1:nCategory]
        block_category_effect_indices = [Int[] for _ in 1:nCategory]

        for marker in eachindex(marker_category_indices)
            categories = marker_category_indices[marker]
            effect_indices = Vector{Int}(undef, length(categories))
            for position in eachindex(categories)
                category = categories[position]
                category_counts[category] += 1
                effect_index = category_counts[category]
                effect_indices[position] = effect_index
                push!(block_category_marker_indices[category], marker)
                push!(block_category_effect_indices[category], effect_index)
            end
            marker_effect_indices[marker] = effect_indices
        end

        block.annotation_effect_indices = marker_effect_indices
        block.category_marker_indices = block_category_marker_indices
        block.category_effect_indices = block_category_effect_indices
    end

    return category_counts
end

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
        block.annotation_category_indices, block.annotation_category_values = build_annotation_category_entries(block.annotation_mask, block.anno_matrix)
    end
    assign_annotation_effect_indices!(blocks, size(first(blocks).annotation_mask, 2))
    return nothing
end

function prepare_marker_probit_tree_block_state!(blocks)
    for block in blocks
        block.xpx, block.x_arrays = build_block_designs(block, 0)
        block.annotation_mask = build_all_marker_mask(block)
        block.annotation_category_indices, block.annotation_category_values = build_annotation_category_entries(block.annotation_mask, Float64.(block.annotation_mask))
    end
    assign_annotation_effect_indices!(blocks, 1)
    annotation_design = build_global_annotation_design(blocks; add_intercept=true)
    return annotation_design
end