using Test
using SBayesAPP

@testset "Non-MPI smoke run" begin
    root = SBayesAPP.repo_root()
    mktempdir() do analysis_dir
        config = SBayesAPP.NonMPIConfig(
            string(joinpath(root, "example", "SBayesAPP_input_first10blks"), "/"),
            string(analysis_dir, "/"),
            2,
            42,
            1,
            "annotation_df.txt",
            "anno_matrix_dict",
            1,
            "XXX",
            "XXX",
            string(joinpath(root, "example", "ST_res"), "/"),
            300000,
            300000,
            1,
            false,
        )

        mkpath(config.analysis_path)
        SBayesAPP.run_nonmpi(config)

        @test isfile(joinpath(analysis_dir, "annotationName.txt"))
        @test isfile(joinpath(analysis_dir, "MCMC_samples_pi.txt"))
        @test isfile(joinpath(analysis_dir, "estPi1.txt"))
        @test isfile(joinpath(analysis_dir, "estR.txt"))
    end
end 