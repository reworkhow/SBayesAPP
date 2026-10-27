@testset "Two-trait marker-probit-tree BayesR mathematics" begin
    design = [1.0 0.0; 1.0 1.0; 1.0 -1.0]
    start_pi = Dict(
        (0.0, 0.0) => 0.7,
        (1.0, 1.0) => 0.1,
        (1.0, 0.0) => 0.1,
        (0.0, 1.0) => 0.1,
    )
    state = SBayesAPP.initialize_marker_probit_tree_state(
        design,
        start_pi;
        prior_model=:marker_probit_tree_bayesr,
    )

    @test size(state.coefficients) == (2, 7)
    @test size(state.snp_pi) == (3, 16)
    @test all(state.snp_pi .== 1 / 16)
    @test length(unique(
        SBayesAPP.marker_probit_tree_bayesr_state_index(pair)
        for pair in SBayesAPP.marker_probit_tree_bayesr_state_order()
    )) == 16

    classes = [
        Float64[1, 2, 1, 2, 3, 4],
        Float64[1, 1, 2, 2, 3, 4],
    ]
    responses, active = SBayesAPP.marker_probit_tree_bayesr_step_indicators(classes)
    @test responses[1] == [0, 1, 1, 1, 1, 1]
    @test responses[2] == [0, 0, 0, 1, 1, 1]
    @test responses[3] == [0, 1, 0, 0, 0, 0]
    @test responses[4] == [0, 0, 0, 0, 1, 1]
    @test responses[5] == [0, 0, 0, 0, 0, 1]
    @test responses[6] == [0, 0, 0, 0, 1, 1]
    @test responses[7] == [0, 0, 0, 0, 0, 1]
    @test active[2] == [2, 3, 4, 5, 6]
    @test active[3] == [2, 3]
    @test active[4] == [2, 4, 5, 6]
    @test active[6] == [3, 4, 5, 6]

    fill!(state.mu, 0.0)
    SBayesAPP.rebuild_marker_probit_tree_bayesr_priors!(state)
    @test vec(sum(state.snp_pi, dims=2)) ≈ ones(3)
    @test state.snp_pi[1, SBayesAPP.marker_probit_tree_bayesr_state_index((1, 1))] ≈ 0.5
    @test state.snp_pi[1, SBayesAPP.marker_probit_tree_bayesr_state_index((2, 1))] ≈ 0.0625
    @test state.snp_pi[1, SBayesAPP.marker_probit_tree_bayesr_state_index((1, 2))] ≈ 0.0625
    @test state.snp_pi[1, SBayesAPP.marker_probit_tree_bayesr_state_index((2, 2))] ≈ 0.0625
    summary = SBayesAPP.marker_probit_tree_summary_dict(state)
    @test summary[(0.0, 0.0)] ≈ 0.5
    @test summary[(1.0, 0.0)] ≈ 0.125
    @test summary[(0.0, 1.0)] ≈ 0.125
    @test summary[(1.0, 1.0)] ≈ 0.25

    log_weights = zeros(4)
    means = zeros(4)
    variances = zeros(4)
    Ginv = [2.0 0.25; 0.25 3.0]
    Rinv = [4.0 0.5; 0.5 5.0]
    SBayesAPP.bayesr_conditional_class_parameters!(
        log_weights,
        means,
        variances,
        state,
        1,
        1,
        3,
        0.4,
        1.2,
        2.5,
        Ginv,
        Rinv,
    )
    @test variances[1] ≈ 1 / Ginv[1, 1]
    @test means[1] ≈ -Ginv[1, 2] * 0.4 / Ginv[1, 1]
    probabilities = zeros(4)
    sampled_class = SBayesAPP.sample_log_categorical!(probabilities, log_weights)
    expected_probabilities = exp.(log_weights .- maximum(log_weights))
    expected_probabilities ./= sum(expected_probabilities)
    @test sampled_class in 1:4
    @test sum(probabilities) ≈ 1.0
    @test probabilities ≈ expected_probabilities
    @test SBayesAPP.marker_probit_tree_bayesr_scale(1) == 0.0
    @test SBayesAPP.marker_probit_tree_bayesr_scale(4) == 1.0
end
