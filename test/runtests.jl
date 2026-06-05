using Test
using SBayesAPP
using CSV
using DataFrames
using DelimitedFiles
using JLD2

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
            "XXX",
            "XXX",
            st_dir,
            1,
            300000,
            300000,
            false,
            burnin=0,
            report_pleiotropic_qtl_effect_matrix=false,
            output_mcmc_delta=false,
        )

        SBayesAPP.run_nonmpi(config)

        @test isfile(joinpath(analysis_dir, "annotationName.txt"))
        @test isfile(joinpath(analysis_dir, "MCMC_samples_pi.txt"))
        @test isfile(joinpath(analysis_dir, "estPi1.txt"))
        @test isfile(joinpath(analysis_dir, "estR.txt"))
        @test isfile(joinpath(analysis_dir, "estGcor.txt"))
        @test isfile(joinpath(analysis_dir, "estGcor_total.txt"))
        @test isfile(joinpath(analysis_dir, "last_sample_delta", "last_sample_delta1_rank0.txt"))
        @test !isfile(joinpath(analysis_dir, "MCMC_samples_marker_effects_variance.txt"))
        @test !isfile(joinpath(analysis_dir, "mcmcAtruecor_c.txt"))
        @test !isfile(joinpath(analysis_dir, "mcmcGcor_c.txt"))
        @test !isfile(joinpath(analysis_dir, "mcmcGcov_c.txt"))
        @test !isfile(joinpath(analysis_dir, "mcmcGcor_total.txt"))
        @test !isfile(joinpath(analysis_dir, "mcmcGcov_total.txt"))
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
            "XXX",
            "XXX",
            st_dir,
            1,
            300000,
            300000,
            false,
            burnin=0,
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
        @test !isfile(joinpath(analysis_dir, "mcmcGcor_total.txt"))
        @test !isfile(joinpath(analysis_dir, "mcmcGcov_total.txt"))
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

@testset "Build non-MPI input dicts" begin
    mktempdir() do tmpdir
        ldroot = joinpath(tmpdir, "ld")
        mkpath(joinpath(ldroot, "readableFiles"))
        outdir = joinpath(tmpdir, "out")

        ldinfo = DataFrame(
            Chrom=[1, 1, 1],
            ID=["rs1", "rs2", "rs3"],
            Index=[1, 2, 3],
            GenPos=[0.0, 0.0, 0.0],
            PhysPos=[10, 20, 30],
            A1=["A", "G", "T"],
            A2=["C", "A", "C"],
            A1Freq=[0.1, 0.2, 0.3],
            N=[100, 100, 100],
            Block=[1, 1, 2],
        )
        CSV.write(joinpath(ldroot, "snp.info"), ldinfo)

        writedlm(joinpath(ldroot, "readableFiles", "block1.lambda.csv"), [4.0, 9.0])
        writedlm(joinpath(ldroot, "readableFiles", "block2.lambda.csv"), [16.0])
        writedlm(joinpath(ldroot, "readableFiles", "block1.U.csv"), [1.0 0.0; 0.0 1.0], ',')
        writedlm(joinpath(ldroot, "readableFiles", "block2.U.csv"), [1.0], ',')

        trait1 = DataFrame(
            SNP=["rs1", "rs2", "rs3"],
            A1=["A", "A", "T"],
            A2=["C", "G", "C"],
            freq=[0.1, 0.8, 0.3],
            b=[0.2, 0.3, 0.4],
            se=[0.1, 0.2, 0.1],
            p=[0.01, 0.02, 0.03],
            N=[100.0, 120.0, 130.0],
        )
        trait2 = DataFrame(
            SNP=["rs1", "rs2", "rs3"],
            A1=["A", "A", "T"],
            A2=["C", "G", "C"],
            freq=[0.1, 0.8, 0.3],
            b=[0.1, 0.2, 0.5],
            se=[0.1, 0.2, 0.1],
            p=[0.04, 0.05, 0.06],
            N=[90.0, 110.0, 140.0],
        )
        annot = DataFrame(SNP=["rs1", "rs2", "rs3"], annot1=[1, 0, 0], annot2=[0, 1, 0], Rest=[0, 0, 1])

        trait1_path = joinpath(tmpdir, "trait1.ma")
        trait2_path = joinpath(tmpdir, "trait2.ma")
        annot_path = joinpath(tmpdir, "annotation_df.txt")
        CSV.write(trait1_path, trait1)
        CSV.write(trait2_path, trait2)
        CSV.write(annot_path, annot)

        result = SBayesAPP.build_nonmpi_input_dicts(outdir, ldroot, trait1_path, trait2_path, "snp.info", annot_path)
        block_data = SBayesAPP.load_nonmpi_block_data(string(dirname(result.transformed_x_file), "/"), result.annot_dict)

        @test block_data.nblk == 2
        @test block_data.nsnp == 3
        @test isfile(result.annotation_file)
        @test isfile(result.transformed_x_file)
        @test isfile(result.transformed_y_file)
        @test isfile(result.blk_snps_index_file)
        @test isfile(result.nGWAS_file)
        @test isfile(result.annotation_dict_file)
        @test block_data.blocks[1].n_gwas == [110.0, 100.0]

        run_annot_file = joinpath(dirname(result.transformed_x_file), "annotation_df.txt")
        cp(result.annotation_file, run_annot_file; force=true)
        st_dir = string(joinpath(SBayesAPP.repo_root(), "example", "ST_res"), "/")
        dense_dir = joinpath(tmpdir, "dense_output")
        sparse_dir = joinpath(tmpdir, "sparse_output")

        function validation_config(analysis_dir)
            return SBayesAPP.NonMPIConfig(
                string(dirname(result.transformed_x_file), "/"),
                string(analysis_dir, "/"),
                4,
                31415,
                "annotation_df.txt",
                result.annot_dict,
                "XXX",
                "XXX",
                st_dir,
                1,
                100,
                100,
                false;
                burnin=0,
                estimate_vare=false,
                estimate_vara=false,
                estGscale_iter=2,
                report_pleiotropic_qtl_effect_matrix=false,
                output_mcmc_delta=false,
            )
        end

        SBayesAPP.run_nonmpi(validation_config(dense_dir))
        SBayesAPP.run_nonmpi(validation_config(sparse_dir))

        comparable_outputs = [
            "MCMC_samples_pi.txt",
            "MCMC_samples_beta_effects_variance.txt",
            "MCMC_samples_genetic_effects_variance.txt",
            "MCMC_samples_total_genetic_effects_variance.txt",
            "estGtotal.txt",
            joinpath("last_sample_delta", "last_sample_delta1_rank0.txt"),
            "last_mcmc_betaArray1.rank0.txt",
            "last_mcmc_betaArray2.rank0.txt",
            "meanAlpha1.rank0.txt",
            "meanAlpha2.rank0.txt",
            joinpath("last_sample_delta", "last_sample_delta2_rank0.txt"),
        ]
        for relative_path in comparable_outputs
            @test read(joinpath(sparse_dir, relative_path), String) == read(joinpath(dense_dir, relative_path), String)
        end
    end
end