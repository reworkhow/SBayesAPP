using Distributions: Chisq, Normal, cdf, truncated
using Random: randn

const MARKER_PROBIT_TREE_STEP_LABELS = (
    "step1_zero_vs_active",
    "step2_11_vs_singleton",
    "step3_10_vs_01",
)

marker_probit_tree_step_labels() = collect(MARKER_PROBIT_TREE_STEP_LABELS)

function marker_probit_tree_state_index(state)
    key = pi_key(state)
    for (index, state_key) in enumerate(pi_key_order())
        if key == state_key
            return index
        end
    end
    error("marker_probit_tree expects binary two-trait states, got $state")
end

function marker_probit_tree_start_row(start_pi::AbstractDict)
    row = zeros(Float64, 4)
    for (state, prob) in start_pi
        row[marker_probit_tree_state_index(state)] = Float64(prob)
    end
    shared_index = marker_probit_tree_state_index((1.0, 1.0))
    trait1_index = marker_probit_tree_state_index((1.0, 0.0))
    trait2_index = marker_probit_tree_state_index((0.0, 1.0))
    trait1_active_mass = row[shared_index] + row[trait1_index]
    trait2_active_mass = row[shared_index] + row[trait2_index]
    shared_mass = row[shared_index]
    isapprox(sum(row), 1.0; atol=1e-8) || error("marker_probit_tree startup Pi must sum to 1.")
    trait1_active_mass > 0.0 || error("marker_probit_tree startup Pi requires positive trait-1 active mass.")
    trait2_active_mass > 0.0 || error("marker_probit_tree startup Pi requires positive trait-2 active mass.")
    shared_mass > 0.0 || error("marker_probit_tree startup Pi requires positive shared-state mass.")
    return row
end

function initialize_marker_probit_tree_state(design_matrix::AbstractMatrix, start_pi::AbstractDict)
    matrix = Matrix{Float64}(design_matrix)
    nmarker, nfeature = size(matrix)
    start_row = marker_probit_tree_start_row(start_pi)
    coefficients = zeros(Float64, nfeature, 3)

    return MarkerProbitTreeState(
        matrix,
        coefficients,
        zeros(Float64, nfeature, 3),
        zeros(Float64, nfeature, 3),
        ones(Float64, 3),
        zeros(Float64, nmarker, 3),
        zeros(Float64, nmarker, 3),
        fill(-Inf, nmarker, 3),
        fill(Inf, nmarker, 3),
        repeat(reshape(start_row, 1, :), nmarker, 1),
    )
end

function marker_probit_tree_step_indicators(deltaArray::AbstractVector)
    d1 = Int.(deltaArray[1])
    d2 = Int.(deltaArray[2])
    states = Vector{Int}(undef, length(d1))

    for idx in eachindex(d1, d2)
        if d1[idx] == 0 && d2[idx] == 0
            states[idx] = 1
        elseif d1[idx] == 1 && d2[idx] == 1
            states[idx] = 2
        elseif d1[idx] == 1 && d2[idx] == 0
            states[idx] = 3
        else
            states[idx] = 4
        end
    end

    z1 = Int.(states .!= 1)
    z2 = Int.(states .== 2)
    z3 = Int.(states .== 3)
    active_sets = (
        collect(eachindex(states)),
        findall(!iszero, z1),
        findall(state -> state == 3 || state == 4, states),
    )
    return (z1, z2, z3), active_sets
end

function sample_marker_probit_tree_liabilities!(
    liability::AbstractVector,
    mu::AbstractVector,
    lower::AbstractVector,
    upper::AbstractVector,
    response::AbstractVector{<:Integer},
)
    for idx in eachindex(response, liability, mu, lower, upper)
        if response[idx] == 0
            lower[idx] = -Inf
            upper[idx] = 0.0
        else
            lower[idx] = 0.0
            upper[idx] = Inf
        end
        liability[idx] = rand(truncated(Normal(mu[idx], 1.0), lower[idx], upper[idx]))
    end
    return nothing
end

function update_marker_probit_tree_coefficients!(
    coefficients::AbstractVector,
    design_matrix::AbstractMatrix,
    latent_residual::AbstractVector,
    coef_prior_var::Real,
)
    nobs = size(design_matrix, 1)
    old_intercept = coefficients[1]

    rhs = sum(latent_residual) + nobs * old_intercept
    inv_lhs = 1.0 / nobs
    ahat = inv_lhs * rhs
    coefficients[1] = randn() * sqrt(inv_lhs) + ahat
    latent_residual .+= old_intercept - coefficients[1]

    if size(design_matrix, 2) > 1
        for feature in 2:size(design_matrix, 2)
            old_value = coefficients[feature]
            x_feature = view(design_matrix, :, feature)
            anno_diag = dot(x_feature, x_feature)
            inv_lhs = 1.0 / (anno_diag + 1.0 / coef_prior_var)
            ahat = inv_lhs * (dot(x_feature, latent_residual) + anno_diag * old_value)
            coefficients[feature] = randn() * sqrt(inv_lhs) + ahat
            latent_residual .+= x_feature .* (old_value - coefficients[feature])
        end
    end
    return nothing
end

function sample_marker_probit_tree_step_variance!(state::MarkerProbitTreeState, step::Int)
    coeffs = view(state.coefficients, :, step)
    n_random_coef = length(coeffs) - 1
    state.variance[step] = (sum(abs2, view(coeffs, 2:length(coeffs))) + 2.0) /
                           rand(Chisq(n_random_coef + 2.0))
    return nothing
end

function sample_marker_probit_tree_step!(
    state::MarkerProbitTreeState,
    step::Int,
    response::AbstractVector{<:Integer},
    active::AbstractVector{<:Integer},
)
    isempty(active) && return nothing

    coeffs = view(state.coefficients, :, step)
    state.mu[:, step] .= state.design_matrix * coeffs
    state.lower_bound[:, step] .= -Inf
    state.upper_bound[:, step] .= Inf

    active_design = state.design_matrix[active, :]
    mu_active = view(state.mu, active, step)
    liability_active = view(state.liability, active, step)
    lower_active = view(state.lower_bound, active, step)
    upper_active = view(state.upper_bound, active, step)
    response_active = view(response, active)

    sample_marker_probit_tree_liabilities!(
        liability_active,
        mu_active,
        lower_active,
        upper_active,
        response_active,
    )

    latent_residual = liability_active .- mu_active
    update_marker_probit_tree_coefficients!(coeffs, active_design, latent_residual, state.variance[step])

    if size(active_design, 2) > 1
        sample_marker_probit_tree_step_variance!(state, step)
    end

    state.mu[:, step] .= state.design_matrix * coeffs
    return nothing
end

function rebuild_marker_probit_tree_priors!(state::MarkerProbitTreeState)
    probs = clamp.(cdf.(Normal(), state.mu), eps(Float64), 1 - eps(Float64))
    p1 = probs[:, 1]
    p2 = probs[:, 2]
    p3 = probs[:, 3]
    state.snp_pi[:, 1] .= 1 .- p1
    state.snp_pi[:, 2] .= p1 .* p2
    state.snp_pi[:, 3] .= p1 .* (1 .- p2) .* p3
    state.snp_pi[:, 4] .= p1 .* (1 .- p2) .* (1 .- p3)
    return nothing
end

function marker_probit_tree_summary_dict(state::MarkerProbitTreeState)
    mean_row = vec(mean(state.snp_pi, dims=1))
    return Dict(state_key => mean_row[index] for (index, state_key) in enumerate(pi_key_order()))
end

function record_marker_probit_tree_coefficient_moments!(state::MarkerProbitTreeState, iIter::Real)
    state.mean_coefficients .+= (state.coefficients .- state.mean_coefficients) .* iIter
    state.mean_coefficients2 .+= ((state.coefficients .^ 2) .- state.mean_coefficients2) .* iIter
    return nothing
end

function update_marker_probit_tree_priors!(state::MarkerProbitTreeState, deltaArray::AbstractVector)
    responses, active_sets = marker_probit_tree_step_indicators(deltaArray)
    for step in 1:3
        sample_marker_probit_tree_step!(state, step, responses[step], active_sets[step])
    end
    rebuild_marker_probit_tree_priors!(state)
    return marker_probit_tree_summary_dict(state)
end

function log_marker_probit_tree_state_prior(state::MarkerProbitTreeState, marker_index::Int, marker_state)
    state_index = marker_probit_tree_state_index(marker_state)
    return log(state.snp_pi[marker_index, state_index])
end

function log_marker_state_prior(
    annotation_prior_model::Symbol,
    Pi,
    marker_probit_tree_state,
    category::Int,
    marker_index::Int,
    marker_state,
)
    if annotation_prior_model == :marker_probit_tree
        marker_probit_tree_state === nothing && error("marker_probit_tree state is required for marker-specific priors.")
        return log_marker_probit_tree_state_prior(marker_probit_tree_state, marker_index, marker_state)
    end
    return log(Pi[category][pi_key(marker_state)])
end