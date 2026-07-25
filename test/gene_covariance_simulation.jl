using Test
using Random
using LinearAlgebra
using DelimitedFiles
using JLD2
using CSV
using DataFrames
using SBayesAPP

function write_gene_recovery_simulation(root; sample_size=250_000, seed=90210)
    data_dir = joinpath(root, "data")
    st_dir = joinpath(root, "st")
    mkpath(data_dir)
    mkpath(joinpath(st_dir, "Trait1"))
    mkpath(joinpath(st_dir, "Trait2"))

    snp_ids = Dict(
        1 => ["snp1", "snp2", "snp3"],
        2 => ["snp4", "snp5"],
    )
    transformed_x = Dict(
        1 => [1.0 0.25 0.0; 0.15 1.0 0.20; 0.0 0.10 0.90],
        2 => [1.0 0.20; 0.30 1.10],
    )
    alpha_trait1 = Dict(1 => [0.35, -0.25, 0.22], 2 => [0.30, -0.28])
    alpha_trait2 = Dict(1 => [0.25, 0.20, -0.18], 2 => [0.24, 0.19])

    rng = MersenneTwister(seed)
    residual_sd = inv(sqrt(sample_size))
    transformed_y = Dict(
        block_id => [
            transformed_x[block_id] * alpha_trait1[block_id] +
            residual_sd * randn(rng, size(transformed_x[block_id], 1)),
            transformed_x[block_id] * alpha_trait2[block_id] +
            residual_sd * randn(rng, size(transformed_x[block_id], 1)),
        ]
        for block_id in keys(transformed_x)
    )

    jldsave(
        joinpath(data_dir, "TransformedX_dict.jld2");
        my_TransformedX_dict=transformed_x,
    )
    jldsave(
        joinpath(data_dir, "TransformedY_dict.jld2");
        my_TransformedY_dict=transformed_y,
    )
    jldsave(
        joinpath(data_dir, "blkSNPsIndex_dict.jld2");
        my_blkSNPsIndex_dict=Dict(1 => collect(1:3), 2 => collect(1:2)),
    )
    jldsave(
        joinpath(data_dir, "nGWAS_dict.jld2");
        my_nGWAS_dict=Dict(
            1 => [Float64(sample_size), Float64(sample_size)],
            2 => [Float64(sample_size), Float64(sample_size)],
        ),
    )
    jldsave(
        joinpath(data_dir, "annotation.jld2");
        my_anno_matrix_dict=Dict(1 => ones(3, 1), 2 => ones(2, 1)),
        my_snp_ids_dict=snp_ids,
    )
    writedlm(joinpath(data_dir, "blkIDs.txt"), [1, 2], ',')

    all_snps = vcat(snp_ids[1], snp_ids[2])
    CSV.write(
        joinpath(data_dir, "annotation_df.txt"),
        DataFrame(SNP=all_snps, category_x=ones(Int, length(all_snps))),
    )
    mapping_path = joinpath(data_dir, "snp_gene_map.tsv")
    open(mapping_path, "w") do io
        println(io, "SNP\tGENE")
        println(io, "snp1\tgeneA")
        println(io, "snp1\tgeneShared")
        println(io, "snp2\tgeneA")
        println(io, "snp3\tgeneB")
        println(io, "snp4\tgeneA")
        println(io, "snp5\tgeneB")
    end

    writedlm(joinpath(st_dir, "Trait1", "mean_varg_total.txt"), [0.5])
    writedlm(joinpath(st_dir, "Trait2", "mean_varg_total.txt"), [0.5])

    return (
        data_dir=data_dir,
        st_dir=st_dir,
        sample_size=sample_size,
        transformed_x=transformed_x,
        alpha_trait1=alpha_trait1,
        alpha_trait2=alpha_trait2,
        genes=Dict(
            "geneA" => Dict(1 => [1, 2], 2 => [1]),
            "geneB" => Dict(1 => [3], 2 => [2]),
            "geneShared" => Dict(1 => [1]),
        ),
        multiplicity=Dict(
            "snp1" => 2,
            "snp2" => 1,
            "snp3" => 1,
            "snp4" => 1,
            "snp5" => 1,
        ),
        snp_ids=snp_ids,
    )
end

function simulated_gene_covariance(simulation, mode::Symbol)
    truth = Dict{String,Float64}()
    for (gene, block_markers) in simulation.genes
        covariance = 0.0
        for (block_id, marker_indices) in block_markers
            weights = [
                mode == :fractional ?
                inv(sqrt(simulation.multiplicity[simulation.snp_ids[block_id][marker]])) :
                1.0
                for marker in marker_indices
            ]
            predictor1 = simulation.transformed_x[block_id][:, marker_indices] *
                         (weights .* simulation.alpha_trait1[block_id][marker_indices])
            predictor2 = simulation.transformed_x[block_id][:, marker_indices] *
                         (weights .* simulation.alpha_trait2[block_id][marker_indices])
            covariance += dot(predictor1, predictor2)
        end
        truth[gene] = covariance
    end
    return truth
end

function run_gene_recovery_simulation(simulation, output_dir, mode::Symbol)
    config = SBayesAPP.NonMPIConfig(
        string(simulation.data_dir, "/"),
        string(output_dir, "/"),
        600,
        8675309,
        "annotation_df.txt",
        "annotation",
        "unused",
        "unused",
        string(simulation.st_dir, "/"),
        2,
        simulation.sample_size,
        simulation.sample_size,
        false;
        burnin=200,
        annotation_prior_model=:group_dirichlet,
        estimate_vare=false,
        estimate_vara=false,
        estimate_pi=false,
        estimate_Gscale=false,
        report_pleiotropic_qtl_effect_matrix=false,
        output_mcmc_delta=false,
        snp_gene_map="snp_gene_map.tsv",
        gene_category="category_x",
        gene_mapping_weight=mode,
    )
    redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            SBayesAPP.run_nonmpi(config)
        end
    end
    return CSV.read(
        joinpath(output_dir, "gene_genetic_covariance_category_x_$(mode).txt"),
        DataFrame;
        delim='\t',
    )
end

@testset "Gene covariance posterior recovery simulation" begin
    mktempdir() do tmpdir
        simulation = write_gene_recovery_simulation(tmpdir)
        for mode in (:full, :fractional)
            output = run_gene_recovery_simulation(
                simulation,
                joinpath(tmpdir, String(mode)),
                mode,
            )
            truth = simulated_gene_covariance(simulation, mode)
            estimated = Dict(
                String(row.GENE) => row.POSTERIOR_MEAN_COV
                for row in eachrow(output)
            )

            @test output.GENE == ["geneA", "geneB", "geneShared"]
            @test output.EFFECTIVE_SNP_COUNT ==
                  (mode == :full ? [3.0, 2.0, 1.0] : [2.5, 2.0, 0.5])
            @test all(isfinite, output.POSTERIOR_MEAN_COV)
            @test all(isfinite, output.POSTERIOR_SD_COV)
            absolute_errors = [
                abs(estimated[gene] - truth[gene])
                for gene in output.GENE
            ]
            @info "Gene covariance recovery" mode genes=Tuple(output.GENE) truth=Tuple(
                round(truth[gene]; digits=6) for gene in output.GENE
            ) estimated=Tuple(
                round(value; digits=6) for value in output.POSTERIOR_MEAN_COV
            ) absolute_errors=Tuple(
                round(value; digits=6) for value in absolute_errors
            )
            @test maximum(absolute_errors) < 0.01
        end
    end
end
