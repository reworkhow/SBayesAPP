using Test
using SBayesAPP

@testset "Pi dict IO is consistent" begin
    pi_dict = Dict(
        (0.0, 0.0) => 0.1,
        (1.0, 1.0) => 0.2,
        (1.0, 0.0) => 0.3,
        (0.0, 1.0) => 0.4,
    )

    mktempdir() do tmpdir
        canonical_path = joinpath(tmpdir, "pi_canonical.txt")
        SBayesAPP.write_pi_dict(canonical_path, pi_dict)
        canonical_read = SBayesAPP.read_to_dict(canonical_path)
        for key in ((0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0))
            @test canonical_read[key] == pi_dict[key]
        end

        legacy_tab_path = joinpath(tmpdir, "pi_legacy_tab.txt")
        write(
            legacy_tab_path,
            "[0.0, 0.0]\t0.1\n[1.0, 1.0]\t0.2\n[1.0, 0.0]\t0.3\n[0.0, 1.0]\t0.4\n",
        )
        legacy_tab_read = SBayesAPP.read_to_dict(legacy_tab_path)
        for key in ((0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0))
            @test legacy_tab_read[key] == pi_dict[key]
        end

        legacy_csv_path = joinpath(tmpdir, "pi_legacy_csv.txt")
        write(
            legacy_csv_path,
            "[0.0,0.0],0.1\n[1.0,1.0],0.2\n[1.0,0.0],0.3\n[0.0,1.0],0.4\n",
        )
        legacy_csv_read = SBayesAPP.read_to_dict(legacy_csv_path)
        for key in ((0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0))
            @test legacy_csv_read[key] == pi_dict[key]
        end
    end
end

@testset "Named CLI args parse correctly" begin
    nonmpi_config = SBayesAPP.parse_nonmpi_args([
        "--data_path", "data/",
        "--analysis_path", "out/",
        "--n_iter", "10",
        "--seed", "1",
        "--nrank", "1",
        "--annot_file", "annot.txt",
        "--annot_dict", "annot_dict",
        "--out_freq", "5",
        "--starting_value_dir", "start",
        "--secondary_starting_value_dir", "second",
        "--st_path", "st/",
        "--thin", "2",
        "--n1", "100",
        "--n2", "200",
        "--n_con", "3",
        "--is_continue", "true",
        "--estimate_vare", "false",
        "--estimate_vara", "true",
        "--estimate_pi", "false",
        "--estimate_gscale", "false",
        "--estgscale_iter", "77",
    ])

    @test nonmpi_config.data_path == "data/"
    @test nonmpi_config.analysis_path == "out/"
    @test nonmpi_config.n_con == 3
    @test nonmpi_config.is_continue
    @test !nonmpi_config.estimate_vare
    @test nonmpi_config.estimate_vara
    @test !nonmpi_config.estimate_pi
    @test !nonmpi_config.estimate_Gscale
    @test nonmpi_config.estGscale_iter == 77
    @test occursin("--analysis_path", string(SBayesAPP.build_nonmpi_cmd(nonmpi_config)))

    mpi_config = SBayesAPP.parse_mpi_args([
        "--data_path=data/",
        "--analysis_path=out/",
        "--n_iter=10",
        "--seed=1",
        "--nrank=2",
        "--annot_file=annot.txt",
        "--annot_dict=annot_dict",
        "--out_freq=5",
        "--starting_value_dir=start",
        "--secondary_starting_value_dir=second",
        "--st_path=st/",
        "--thin=2",
        "--n1=100",
        "--n2=200",
        "--n_con=1",
        "--estimate_pi=true",
        "--fixed_hyperparameters=false",
        "--is_continue=true",
        "--chr=22",
    ])

    @test mpi_config.data_path == "data/"
    @test mpi_config.analysis_path == "out/"
    @test mpi_config.n_con == 1
    @test mpi_config.estimate_pi
    @test !mpi_config.fixed_hyperparameters
    @test mpi_config.is_continue
    @test mpi_config.chr == "22"
    @test occursin("--analysis_path", string(SBayesAPP.build_mpi_cmd(mpi_config)))

    @test_throws ErrorException SBayesAPP.parse_nonmpi_args([
        "data/", "out/", "10", "1", "1", "annot.txt", "annot_dict", "5", "start", "second", "st/", "2", "100", "200", "false",
    ])

    @test_throws ErrorException SBayesAPP.parse_mpi_args([
        "data/", "out/", "10", "1", "2", "annot.txt", "annot_dict", "5", "start", "second", "st/", "2", "100", "200", "true",
    ])
end

@testset "Annotation metadata supports mixed types" begin
    mktempdir() do tmpdir
        annot_path = joinpath(tmpdir, "mixed_annotations.txt")
        write(
            annot_path,
            "SNP\tcontinuous1\tcategory1\nrs1\t0.5\t1\nrs2\t0.0\t0\nrs3\t1.5\t1\n",
        )

        metadata = SBayesAPP.load_annotation_metadata(string(tmpdir, "/"), "mixed_annotations.txt"; nCon=1)

        @test metadata.annotationName == ["continuous1", "category1"]
        @test metadata.nLoci_annot == [2, 2]
        @test metadata.nCon == 1
        @test metadata.nCat == 1
        @test metadata.annotationType == ["continue", "category"]
    end
end

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
            1,
            300000,
            300000,
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