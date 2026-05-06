using Test
using SBayesAPP

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
            1,
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