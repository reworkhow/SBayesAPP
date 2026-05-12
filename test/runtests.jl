using Test
using SBayesAPP

include(joinpath(SBayesAPP.repo_root(), "example", "simulate_marker_probit_tree_dataset.jl"))

@testset "Non-MPI smoke run" begin
    root = SBayesAPP.repo_root()
    mktempdir() do tmpdir
        analysis_dir = joinpath(tmpdir, "nested_output")
        data_dir = string(joinpath(root, "example", "SBayesAPP_input_first10blks"), "/")
        st_dir = string(joinpath(root, "example", "ST_res"), "/")
        config = SBayesAPP.NonMPIConfig(
            data_dir,
            string(analysis_dir, "/"),
            2,
            42,
            "annotation_df.txt",
            "anno_matrix_dict",
            1,
            "XXX",
            "XXX",
            st_dir,
            1,
            300000,
            300000,
            false,
            report_pleiotropic_qtl_effect_matrix=false,
            output_mcmc_delta=false,
        )

        SBayesAPP.run_nonmpi(config)

        @test isfile(joinpath(analysis_dir, "annotationName.txt"))
        @test isfile(joinpath(analysis_dir, "MCMC_samples_pi.txt"))
        @test isfile(joinpath(analysis_dir, "estPi1.txt"))
        @test isfile(joinpath(analysis_dir, "estR.txt"))
        @test isfile(joinpath(analysis_dir, "last_sample_delta", "last_sample_delta1_rank0.txt"))
        @test !isfile(joinpath(analysis_dir, "MCMC_samples_marker_effects_variance.txt"))
        @test !isfile(joinpath(analysis_dir, "mcmcAtruecor_c.txt"))
        @test !isfile(joinpath(analysis_dir, "mcmc_Delta1.rank0.txt"))
        @test !isfile(joinpath(analysis_dir, "estA1.txt"))
    end
end 

@testset "Non-MPI smoke run marker_probit_tree" begin
    root = SBayesAPP.repo_root()
    data_dir = string(joinpath(root, "example", "SBayesAPP_input_first10blks"), "/")
    st_dir = string(joinpath(root, "example", "ST_res"), "/")
    mktempdir() do tmpdir
        analysis_dir = joinpath(tmpdir, "marker_probit_tree_output")
        config = SBayesAPP.NonMPIConfig(
            data_dir,
            string(analysis_dir, "/"),
            2,
            42,
            "annotation_df.txt",
            "anno_matrix_dict",
            1,
            "XXX",
            "XXX",
            st_dir,
            1,
            300000,
            300000,
            false,
            annotation_prior_model=:marker_probit_tree,
            report_pleiotropic_qtl_effect_matrix=false,
            output_mcmc_delta=false,
        )

        SBayesAPP.run_nonmpi(config)

        @test isfile(joinpath(analysis_dir, "annotationName.txt"))
        @test isfile(joinpath(analysis_dir, "annotation_probit_coefficients.txt"))
        @test isfile(joinpath(analysis_dir, "annotation_probit_coefficients_std.txt"))
        @test isfile(joinpath(analysis_dir, "estPi1.txt"))
        @test !isfile(joinpath(analysis_dir, "estPi2.txt"))
        @test isfile(joinpath(analysis_dir, "last_sample_delta", "last_sample_delta1_rank0.txt"))
        @test !isfile(joinpath(analysis_dir, "estBcor.txt"))
        @test !isfile(joinpath(analysis_dir, "estBcor_std.txt"))
        @test !isfile(joinpath(analysis_dir, "estGcor.txt"))
        @test !isfile(joinpath(analysis_dir, "estGcor_std.txt"))
        @test !isfile(joinpath(analysis_dir, "mcmcGcor_c.txt"))
        @test !isfile(joinpath(analysis_dir, "mcmcGcov_c.txt"))
    end
end

@testset "marker_probit_tree synthetic generator" begin
    mktempdir() do tmpdir
        synthetic = simulate_marker_probit_tree_dataset(joinpath(tmpdir, "synthetic_case"))
        metadata = SBayesAPP.load_annotation_metadata(synthetic.data_path, "annotation_df.txt")
        block_data = SBayesAPP.load_nonmpi_block_data(synthetic.data_path, "anno_matrix_dict")

        @test metadata.annotationName == ["annotation1"]
        @test metadata.nCat == 1
        @test metadata.nCon == 0
        @test block_data.nblk == 2
        @test block_data.nsnp == 400
        @test isfile(joinpath(synthetic.truth_path, "annotation_probit_coefficients_truth.txt"))
        @test isfile(joinpath(synthetic.truth_path, "mean_state_probabilities_truth.txt"))
        @test isfile(joinpath(synthetic.truth_path, "realized_state_counts.txt"))
    end
end