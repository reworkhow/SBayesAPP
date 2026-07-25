using Test
using LinearAlgebra
using DelimitedFiles
using JLD2
using CSV
using DataFrames
using SBayesAPP

gene_cli_args() = [
    "--data_path", "data/",
    "--analysis_path", "analysis/",
    "--n_iter", "2",
    "--burnin", "0",
    "--annot_file", "annotation_df.txt",
    "--annot_dict", "anno_matrix_dict",
    "--starting_value_dir", "XXX",
    "--gscale_value_dir", "XXX",
    "--st_path", "st/",
    "--thin", "1",
    "--n1", "100",
    "--n2", "100",
    "--is_continue", "false",
]

function captured_error(f)
    try
        f()
        return ""
    catch error
        return sprint(showerror, error)
    end
end

function gene_test_block(block_id, transformed_x, snp_indices, anno_matrix)
    return SBayesAPP.NonMPIBlockData(
        block_id,
        transformed_x,
        [zeros(size(transformed_x, 1)), zeros(size(transformed_x, 1))],
        snp_indices,
        [100.0, 100.0],
        anno_matrix,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )
end

function write_gene_mapping(path, rows; delim='\t')
    open(path, "w") do io
        println(io, "SNP", delim, "GENE")
        for (snp, gene) in rows
            println(io, snp, delim, gene)
        end
    end
end

@testset "Gene covariance configuration" begin
    config = SBayesAPP.parse_nonmpi_args(gene_cli_args())
    @test (config.snp_gene_map, config.gene_category, config.gene_mapping_weight) ==
          (nothing, nothing, :full)
    @test !occursin("--snp-gene-map", string(SBayesAPP.build_nonmpi_cmd(config)))

    @test_throws ErrorException SBayesAPP.parse_nonmpi_args(
        vcat(gene_cli_args(), ["--snp-gene-map", "map.tsv"]),
    )
    @test_throws ErrorException SBayesAPP.parse_nonmpi_args(
        vcat(gene_cli_args(), ["--gene-category", "category1"]),
    )
    @test_throws ErrorException SBayesAPP.parse_nonmpi_args([
        "--gene-mapping-weight", "fractional",
    ])

    config = SBayesAPP.parse_nonmpi_args(vcat(gene_cli_args(), [
        "--snp-gene-map", "map.tsv",
        "--gene-category", "category1",
        "--gene-mapping-weight", "fractional",
    ]))
    @test (config.snp_gene_map, config.gene_category, config.gene_mapping_weight) ==
          ("map.tsv", "category1", :fractional)

    @test_throws ErrorException SBayesAPP.parse_nonmpi_args(vcat(gene_cli_args(), [
        "--annotation-prior-model", "marker_probit_tree",
        "--snp-gene-map", "map.tsv",
        "--gene-category", "category1",
    ]))
end

function write_gene_runtime_files(tmpdir)
    jldsave(
        joinpath(tmpdir, "TransformedX_dict.jld2");
        my_TransformedX_dict=Dict(1 => [1.0 0.2; 0.2 1.0]),
    )
    jldsave(
        joinpath(tmpdir, "TransformedY_dict.jld2");
        my_TransformedY_dict=Dict(1 => [zeros(2), zeros(2)]),
    )
    jldsave(
        joinpath(tmpdir, "blkSNPsIndex_dict.jld2");
        my_blkSNPsIndex_dict=Dict(1 => [1, 2]),
    )
    jldsave(
        joinpath(tmpdir, "nGWAS_dict.jld2");
        my_nGWAS_dict=Dict(1 => [100.0, 100.0]),
    )
    writedlm(joinpath(tmpdir, "blkIDs.txt"), [1], ',')
end

@testset "Gene covariance SNP-ID compatibility" begin
    mktempdir() do tmpdir
        ldinfo_path = joinpath(tmpdir, "snp.info")
        annotation_path = joinpath(tmpdir, "annotation.tsv")
        CSV.write(ldinfo_path, DataFrame(
            Block=[1, 1],
            ID=["rs_late", "rs_early"],
            Index=[2, 1],
            A1=["A", "C"],
            A2=["T", "G"],
        ); delim='\t')
        CSV.write(annotation_path, DataFrame(
            SNP=["rs_late", "rs_early"],
            category1=[0, 1],
        ); delim='\t')
        write_gene_runtime_files(tmpdir)

        result = SBayesAPP.build_annotation_dict(
            tmpdir,
            annotation_path,
            ldinfo_path;
            output_file_name="new_annotation.jld2",
        )
        annotation_data = JLD2.load(result.output_file)
        @test annotation_data["my_snp_ids_dict"][1] == ["rs_early", "rs_late"]
        @test vec(annotation_data["my_anno_matrix_dict"][1]) == [1.0, 0.0]

        old_path = joinpath(tmpdir, "old_annotation.jld2")
        jldsave(
            old_path;
            my_anno_matrix_dict=annotation_data["my_anno_matrix_dict"],
        )
        old_blocks = SBayesAPP.load_nonmpi_block_data(tmpdir, "old_annotation")
        new_blocks = SBayesAPP.load_nonmpi_block_data(tmpdir, "new_annotation")
        @test old_blocks.nsnp == new_blocks.nsnp == 2
        @test !(:snp_ids in fieldnames(typeof(new_blocks.blocks[1])))
        @test occursin(
            "Rebuild the annotation dictionary",
            captured_error(() -> SBayesAPP.load_gene_snp_ids(old_path)),
        )

        snp_ids = SBayesAPP.load_gene_snp_ids(result.output_file)
        @test isnothing(
            SBayesAPP.validate_gene_snp_id_alignment(new_blocks.blocks, snp_ids),
        )
        snp_ids[1] = ["rs_early"]
        @test_throws ErrorException SBayesAPP.validate_gene_snp_id_alignment(
            new_blocks.blocks,
            snp_ids,
        )
    end
end

function gene_mapping_fixture()
    blocks = [
        gene_test_block(
            1,
            [1.0 0.2 0.0; 0.2 1.0 0.1],
            [1, 2, 3],
            [1.0 0.0; 1.0 0.0; 0.0 1.0],
        ),
        gene_test_block(
            2,
            [1.0 0.3 0.0; 0.3 1.0 0.2],
            [4, 5, 6],
            [1.0 0.0; 0.0 1.0; 1.0 0.0],
        ),
    ]
    SBayesAPP.prepare_block_state!(blocks, 0, ["category", "category"])
    return (
        blocks=blocks,
        snp_ids=Dict(1 => ["A", "B", "C"], 2 => ["D", "E", "F"]),
        metadata=(
            annotationName=["category_x", "other"],
            annotationType=["category", "category"],
        ),
    )
end

@testset "Gene covariance sparse mapping and weighting" begin
    mktempdir() do tmpdir
        fixture = gene_mapping_fixture()
        mapping_rows = [
            ("A", "geneB"),
            ("A", "geneB"),
            ("A", "geneA"),
            ("B", "geneA"),
            ("C", "outsideGene"),
            ("D", "geneA"),
            ("Z", "absentGene"),
        ]
        mapping_path = joinpath(tmpdir, "mapping.tsv")
        write_gene_mapping(mapping_path, mapping_rows)

        full_state = SBayesAPP.initialize_gene_covariance_state(
            fixture.blocks,
            fixture.snp_ids,
            fixture.metadata,
            mapping_path,
            "category_x",
            :full,
        )
        @test full_state.gene_names == ["geneA", "geneB"]
        @test full_state.snp_counts == [3, 1]
        @test full_state.effective_snp_counts == [3.0, 1.0]
        @test length(full_state.per_block_slices[1]) == 2
        @test length(full_state.per_block_slices[2]) == 1

        gene_a_block1 = only(filter(
            slice -> slice.gene_index == 1,
            full_state.per_block_slices[1],
        ))
        gene_a_block2 = only(filter(
            slice -> slice.gene_index == 1,
            full_state.per_block_slices[2],
        ))
        @test gene_a_block1.marker_indices == [1, 2]
        @test gene_a_block1.effect_indices == [1, 2]
        @test gene_a_block2.marker_indices == [1]
        @test gene_a_block2.effect_indices == [3]
        @test full_state.qc.valid_unique_snp_gene_pairs == 4
        @test full_state.qc.valid_category_mapped_snps == 3
        @test full_state.qc.category_snps_without_valid_gene_mapping == 1
        @test full_state.qc.snps_mapped_to_multiple_retained_genes == 1

        fractional_state = SBayesAPP.initialize_gene_covariance_state(
            fixture.blocks,
            fixture.snp_ids,
            fixture.metadata,
            mapping_path,
            "category_x",
            :fractional,
        )
        @test fractional_state.effective_snp_counts == [2.5, 0.5]
        fractional_gene_a_block1 = only(filter(
            slice -> slice.gene_index == 1,
            fractional_state.per_block_slices[1],
        ))
        @test fractional_gene_a_block1.weights ≈ [inv(sqrt(2.0)), 1.0]

        csv_path = joinpath(tmpdir, "mapping.csv")
        write_gene_mapping(csv_path, mapping_rows; delim=',')
        @test SBayesAPP.read_snp_gene_map(csv_path) ==
              SBayesAPP.read_snp_gene_map(mapping_path)

        flexible_path = joinpath(tmpdir, "mapping_with_extra_column.csv")
        open(flexible_path, "w") do io
            println(io, "GENE,SNP,NOTE")
            println(io, "geneA,A,ignored")
            println(io, ",B,missing gene")
        end
        @test SBayesAPP.read_snp_gene_map(flexible_path) == [("A", "geneA")]

        @test_throws ErrorException SBayesAPP.initialize_gene_covariance_state(
            fixture.blocks,
            fixture.snp_ids,
            fixture.metadata,
            mapping_path,
            "missing_category",
            :full,
        )
        continuous_metadata = (
            annotationName=["category_x", "other"],
            annotationType=["continue", "category"],
        )
        @test_throws ErrorException SBayesAPP.initialize_gene_covariance_state(
            fixture.blocks,
            fixture.snp_ids,
            continuous_metadata,
            mapping_path,
            "category_x",
            :full,
        )
        fixture.blocks[1].anno_matrix[1, 1] = 0.5
        @test_throws ErrorException SBayesAPP.initialize_gene_covariance_state(
            fixture.blocks,
            fixture.snp_ids,
            fixture.metadata,
            mapping_path,
            "category_x",
            :full,
        )
    end
end

function covariance_fixture(tmpdir, rows; mapping_weight=:full)
    blocks = [
        gene_test_block(
            1,
            [1.0 0.4; 0.2 1.1],
            [1, 2],
            reshape([1.0, 1.0], 2, 1),
        ),
        gene_test_block(
            2,
            reshape([1.3, 0.5], 2, 1),
            [3],
            reshape([1.0], 1, 1),
        ),
    ]
    SBayesAPP.prepare_block_state!(blocks, 0, ["category"])
    mapping_path = joinpath(tmpdir, "covariance_mapping.tsv")
    write_gene_mapping(mapping_path, rows)
    state = SBayesAPP.initialize_gene_covariance_state(
        blocks,
        Dict(1 => ["A", "B"], 2 => ["D"]),
        (annotationName=["category_x"], annotationType=["category"]),
        mapping_path,
        "category_x",
        mapping_weight,
    )
    return state, blocks
end

@testset "Gene covariance calculation and posterior moments" begin
    mktempdir() do tmpdir
        state, blocks = covariance_fixture(
            tmpdir,
            [("A", "geneA"), ("B", "geneA"), ("D", "geneA")],
        )
        alpha = [
            [[0.7, -0.3, 0.5]],
            [[0.2, 0.6, -0.4]],
        ]
        expected_first =
            dot(
                blocks[1].transformed_x * alpha[1][1][1:2],
                blocks[1].transformed_x * alpha[2][1][1:2],
            ) +
            dot(
                blocks[2].transformed_x * alpha[1][1][3:3],
                blocks[2].transformed_x * alpha[2][1][3:3],
            )
        alpha_before = deepcopy(alpha)
        first_sample = SBayesAPP.accumulate_gene_covariance_sample!(state, blocks, alpha)
        @test first_sample ≈ [expected_first]
        @test alpha == alpha_before

        alpha[1][1] .= [-0.1, 0.8, 0.25]
        alpha[2][1] .= [0.9, -0.2, 0.7]
        expected_second =
            dot(
                blocks[1].transformed_x * alpha[1][1][1:2],
                blocks[1].transformed_x * alpha[2][1][1:2],
            ) +
            dot(
                blocks[2].transformed_x * alpha[1][1][3:3],
                blocks[2].transformed_x * alpha[2][1][3:3],
            )
        second_sample = SBayesAPP.accumulate_gene_covariance_sample!(state, blocks, alpha)
        @test second_sample ≈ [expected_second]
        @test state.posterior_mean_cov ≈ [(expected_first + expected_second) / 2]
        @test state.posterior_second_moment_cov ≈
              [(expected_first^2 + expected_second^2) / 2]
        @test state.retained_sample_count == 2

        fractional_state, fractional_blocks = covariance_fixture(
            tmpdir,
            [("A", "geneA"), ("A", "geneB")];
            mapping_weight=:fractional,
        )
        fractional_alpha = [
            [[0.8, 0.0, 0.0]],
            [[-0.5, 0.0, 0.0]],
        ]
        original_covariance = dot(
            fractional_blocks[1].transformed_x[:, 1] * fractional_alpha[1][1][1],
            fractional_blocks[1].transformed_x[:, 1] * fractional_alpha[2][1][1],
        )
        fractional_sample = SBayesAPP.accumulate_gene_covariance_sample!(
            fractional_state,
            fractional_blocks,
            fractional_alpha,
        )
        @test fractional_sample ≈ fill(original_covariance / 2, 2)
        @test sum(fractional_sample) ≈ original_covariance
    end
end

function gene_output_qc(; retained_genes=2)
    return SBayesAPP.GeneCovarianceQC(
        3,
        2,
        1,
        1,
        1,
        2,
        retained_genes,
    )
end

@testset "Gene covariance output" begin
    mktempdir() do tmpdir
        state = SBayesAPP.GeneCovarianceState(
            ["geneA", "geneB"],
            "red cell / category",
            3,
            :fractional,
            Vector{SBayesAPP.GeneBlockSlice}[],
            [1, 2],
            [0.5, 1.5],
            [-0.2, 0.4],
            [0.08, 0.2],
            4,
            gene_output_qc(),
        )
        paths = SBayesAPP.write_gene_covariance_outputs(tmpdir, state)
        @test basename(paths.main) ==
              "gene_genetic_covariance_red_cell_category_fractional.txt"
        @test basename(paths.qc) ==
              "gene_genetic_covariance_red_cell_category_fractional_qc.txt"

        main_lines = readlines(paths.main)
        @test main_lines[1] ==
              "GENE\tN_CATEGORY_X_SNPS\tEFFECTIVE_SNP_COUNT\tPOSTERIOR_MEAN_COV\tPOSTERIOR_SD_COV"
        @test startswith(main_lines[2], "geneA\t1\t0.5\t-0.2\t")
        @test startswith(main_lines[3], "geneB\t2\t1.5\t0.4\t")
        qc_text = read(paths.qc, String)
        @test occursin("GENE_MAPPING_WEIGHT\tfractional", qc_text)
        @test occursin("RETAINED_MCMC_SAMPLES\t4", qc_text)

        empty_state = SBayesAPP.GeneCovarianceState(
            String[],
            "empty",
            1,
            :full,
            Vector{SBayesAPP.GeneBlockSlice}[],
            Int[],
            Float64[],
            Float64[],
            Float64[],
            0,
            gene_output_qc(retained_genes=0),
        )
        empty_paths = SBayesAPP.write_gene_covariance_outputs(tmpdir, empty_state)
        @test length(readlines(empty_paths.main)) == 1
        @test occursin("RETAINED_GENES\t0", read(empty_paths.qc, String))
    end
end
