using Distributions: Chisq, Normal, cdf, truncated
using Random: randn

const MARKER_PROBIT_TREE_STEP_LABELS = (
    "step1_zero_vs_active",
    "step2_11_vs_singleton",
    "step3_10_vs_01",
)
const MARKER_PROBIT_TREE_BAYESR_STEP_LABELS = (
    "step1_zero_vs_active",
    "step2_11_vs_singleton",
    "step3_10_vs_01",
    "step4_trait1_medium_large_vs_small",
    "step5_trait1_large_vs_medium",
    "step6_trait2_medium_large_vs_small",
    "step7_trait2_large_vs_medium",
)
const MARKER_PROBIT_TREE_BAYESR_GAMMA = (0.0, 0.01, 0.1, 1.0)
const MARKER_PROBIT_TREE_MODELS = (:marker_probit_tree, :marker_probit_tree_bayesr)

marker_probit_tree_step_labels() = collect(MARKER_PROBIT_TREE_STEP_LABELS)
marker_probit_tree_step_labels(state::MarkerProbitTreeState) =
    state.prior_model == :marker_probit_tree_bayesr ?
    collect(MARKER_PROBIT_TREE_BAYESR_STEP_LABELS) :
    marker_probit_tree_step_labels()

is_marker_probit_tree_model(model::Symbol) = model in MARKER_PROBIT_TREE_MODELS
is_marker_probit_tree_bayesr(model::Symbol) = model == :marker_probit_tree_bayesr

function marker_probit_tree_bayesr_state_index(state)
    class1 = Int(state[1])
    class2 = Int(state[2])
    1 <= class1 <= 4 || error("BayesR trait-1 class must be in 1:4, got $class1")
    1 <= class2 <= 4 || error("BayesR trait-2 class must be in 1:4, got $class2")
    return (class1 - 1) * 4 + class2
end

marker_probit_tree_bayesr_state_order() =
    [(class1, class2) for class1 in 1:4 for class2 in 1:4]

marker_probit_tree_bayesr_scale(class::Integer) = sqrt(MARKER_PROBIT_TREE_BAYESR_GAMMA[class])

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

function initialize_marker_probit_tree_state(
    design_matrix::AbstractMatrix,
    start_pi::AbstractDict;
    prior_model::Symbol=:marker_probit_tree,
)
    is_marker_probit_tree_model(prior_model) || error("Unsupported marker-probit-tree model: $prior_model")
    matrix = Matrix{Float64}(design_matrix)
    nmarker, nfeature = size(matrix)
    is_bayesr = is_marker_probit_tree_bayesr(prior_model)
    nsteps = is_bayesr ? 7 : 3
    nstates = is_bayesr ? 16 : 4
    start_row = is_bayesr ? fill(1.0 / nstates, nstates) : marker_probit_tree_start_row(start_pi)
    coefficients = zeros(Float64, nfeature, nsteps)

    return MarkerProbitTreeState(
        prior_model,
        matrix,
        coefficients,
        zeros(Float64, nfeature, nsteps),
        zeros(Float64, nfeature, nsteps),
        ones(Float64, nsteps),
        zeros(Float64, nmarker, nsteps),
        zeros(Float64, nmarker, nsteps),
        fill(-Inf, nmarker, nsteps),
        fill(Inf, nmarker, nsteps),
        repeat(reshape(start_row, 1, :), nmarker, 1),
        zeros(Float64, nstates),
        zeros(Float64, nstates),
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

function marker_probit_tree_bayesr_step_indicators(deltaArray::AbstractVector)
    class1 = Int.(deltaArray[1])
    class2 = Int.(deltaArray[2])
    all(class -> 1 <= class <= 4, class1) || error("BayesR trait-1 classes must be in 1:4.")
    all(class -> 1 <= class <= 4, class2) || error("BayesR trait-2 classes must be in 1:4.")

    active1 = class1 .> 1
    active2 = class2 .> 1
    any_active = active1 .| active2
    shared = active1 .& active2
    singleton = xor.(active1, active2)

    responses = (
        Int.(any_active),
        Int.(shared),
        Int.(active1 .& singleton),
        Int.(class1 .>= 3),
        Int.(class1 .== 4),
        Int.(class2 .>= 3),
        Int.(class2 .== 4),
    )
    active_sets = (
        collect(eachindex(class1)),
        findall(any_active),
        findall(singleton),
        findall(active1),
        findall(class1 .>= 3),
        findall(active2),
        findall(class2 .>= 3),
    )
    return responses, active_sets
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

function rebuild_marker_probit_tree_bayesr_priors!(state::MarkerProbitTreeState)
    # prevents probabilities from becoming exactly 0 or 1
    # l_jh = a_j'\theta_h + epsilon_jh, epsilon_jh ~ N(0,1)
    # Pr(z_jh = 1) = Pr(l_jh >0) = Pr(a_j'\theta_h + epsilon_jh > 0) = Pr(epsilon_jh > -a_j'\theta_h) =\Phi(a_j'\theta_h) 
    probs = clamp.(cdf.(Normal(), state.mu), eps(Float64), 1 - eps(Float64))
    p1, p2, p3 = (view(probs, :, step) for step in 1:3)
    pattern00 = 1 .- p1
    pattern11 = p1 .* p2
    pattern10 = p1 .* (1 .- p2) .* p3
    pattern01 = p1 .* (1 .- p2) .* (1 .- p3)

    q11, q12, q21, q22 = (view(probs, :, step) for step in 4:7)
    #m_t,small = 1 .- qt1
    #m_t,median = qt1 * (1-qt2)
    #m_t,large = qt1 * qt2
    magnitude1 = hcat(1 .- q11, q11 .* (1 .- q12), q11 .* q12)
    magnitude2 = hcat(1 .- q21, q21 .* (1 .- q22), q21 .* q22)

    fill!(state.snp_pi, 0.0)
    state.snp_pi[:, marker_probit_tree_bayesr_state_index((1, 1))] .= pattern00
    for class in 2:4
        state.snp_pi[:, marker_probit_tree_bayesr_state_index((class, 1))] .=
            pattern10 .* view(magnitude1, :, class - 1)
        state.snp_pi[:, marker_probit_tree_bayesr_state_index((1, class))] .=
            pattern01 .* view(magnitude2, :, class - 1)
    end
    for class1 in 2:4, class2 in 2:4
        state.snp_pi[:, marker_probit_tree_bayesr_state_index((class1, class2))] .=
            pattern11 .* view(magnitude1, :, class1 - 1) .* view(magnitude2, :, class2 - 1)
    end
    state.snp_pi ./= sum(state.snp_pi, dims=2)
    return nothing
end

function marker_probit_tree_summary_dict(state::MarkerProbitTreeState)
    if !is_marker_probit_tree_bayesr(state.prior_model)
        mean_row = vec(mean(state.snp_pi, dims=1))
        return Dict(state_key => mean_row[index] for (index, state_key) in enumerate(pi_key_order()))
    end

    mean_row = vec(mean(state.snp_pi, dims=1))
    p00 = mean_row[marker_probit_tree_bayesr_state_index((1, 1))]
    p10 = sum(mean_row[marker_probit_tree_bayesr_state_index((class1, 1))] for class1 in 2:4)
    p01 = sum(mean_row[marker_probit_tree_bayesr_state_index((1, class2))] for class2 in 2:4)
    p11 = sum(
        mean_row[marker_probit_tree_bayesr_state_index((class1, class2))]
        for class1 in 2:4 for class2 in 2:4
    )
    return Dict((0.0, 0.0) => p00, (1.0, 1.0) => p11, (1.0, 0.0) => p10, (0.0, 1.0) => p01)
end

function record_marker_probit_tree_coefficient_moments!(state::MarkerProbitTreeState, iIter::Real)
    state.mean_coefficients .+= (state.coefficients .- state.mean_coefficients) .* iIter
    state.mean_coefficients2 .+= ((state.coefficients .^ 2) .- state.mean_coefficients2) .* iIter
    mean_snp_pi = vec(mean(state.snp_pi, dims=1))
    state.mean_snp_pi .+= (mean_snp_pi .- state.mean_snp_pi) .* iIter
    state.mean_snp_pi2 .+= ((mean_snp_pi .^ 2) .- state.mean_snp_pi2) .* iIter
    return nothing
end

function update_marker_probit_tree_priors!(state::MarkerProbitTreeState, deltaArray::AbstractVector)
    is_bayesr = is_marker_probit_tree_bayesr(state.prior_model)
    # 0/1 responses, row index for active sets
    responses, active_sets = is_bayesr ?
        marker_probit_tree_bayesr_step_indicators(deltaArray) :
        marker_probit_tree_step_indicators(deltaArray)
    for step in eachindex(responses)
        sample_marker_probit_tree_step!(state, step, responses[step], active_sets[step])
    end
    if is_bayesr
        rebuild_marker_probit_tree_bayesr_priors!(state)
    else
        rebuild_marker_probit_tree_priors!(state)
    end
    return marker_probit_tree_summary_dict(state)
end

function log_marker_probit_tree_state_prior(state::MarkerProbitTreeState, marker_index::Int, marker_state)
    state_index = is_marker_probit_tree_bayesr(state.prior_model) ?
        marker_probit_tree_bayesr_state_index(marker_state) : # get state_index from the 16 states
        marker_probit_tree_state_index(marker_state)          # get state_index from the 4 states
    return log(state.snp_pi[marker_index, state_index])
end

function bayesr_conditional_class_parameters!(
    log_weights::AbstractVector,
    means::AbstractVector,
    variances::AbstractVector,
    state::MarkerProbitTreeState,
    marker_index::Int,
    trait::Int,
    other_class::Int,
    other_beta::Real,
    projected_residual::Real,
    xpx::Real,
    Ginv::AbstractMatrix,
    Rinv::AbstractMatrix,
)
    other = trait == 1 ? 2 : 1
    other_scale = marker_probit_tree_bayesr_scale(other_class)
    # when the class for trait = candidate_class
    for candidate_class in 1:4
        scale = marker_probit_tree_bayesr_scale(candidate_class)
        precision = Ginv[trait, trait] + xpx * scale^2 * Rinv[trait, trait] # Cj,11
        cross_precision = Ginv[trait, other] + xpx * scale * other_scale * Rinv[trait, other] # Cj,12
        rhs = scale * projected_residual - cross_precision * other_beta
        variance = 1.0 / precision
        mean_effect = rhs * variance

        joint_state = trait == 1 ? (candidate_class, other_class) : (other_class, candidate_class)
        log_weights[candidate_class] =
            log_marker_probit_tree_state_prior(state, marker_index, joint_state) -
            0.5 * log(precision) + 0.5 * rhs * mean_effect

        means[candidate_class] = mean_effect
        variances[candidate_class] = variance
    end
    return nothing
end

function sample_log_categorical!(probabilities::AbstractVector, log_weights::AbstractVector)
    max_log_weight = maximum(log_weights)
    total = 0.0
    for class in eachindex(log_weights, probabilities)
        # prevents numerical overflow for exp
        # exp(l_a - l_max)/exp(l_b - l_max) = exp(l_a)/exp(l_b)
        #prob_i = exp(l_i)/sum_j exp(l_j)
        probabilities[class] = exp(log_weights[class] - max_log_weight) 
        total += probabilities[class] 
    end
    isfinite(total) && total > 0.0 || error("Invalid BayesR class probabilities.")
    for class in eachindex(probabilities)
        probabilities[class] /= total
    end
    draw = rand()
    cumulative = 0.0
    for class in eachindex(probabilities)
        cumulative += probabilities[class]
        draw <= cumulative && return class
    end
    return lastindex(probabilities)
end

function log_marker_state_prior(
    annotation_prior_model::Symbol,
    Pi,
    marker_probit_tree_state,
    category::Int,
    marker_index::Int,
    marker_state,
)
    if is_marker_probit_tree_model(annotation_prior_model)
        marker_probit_tree_state === nothing && error("marker_probit_tree state is required for marker-specific priors.")
        return log_marker_probit_tree_state_prior(marker_probit_tree_state, marker_index, marker_state)
    end
    return log(Pi[category][pi_key(marker_state)])
end
