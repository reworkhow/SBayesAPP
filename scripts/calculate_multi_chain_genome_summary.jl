using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using SBayesAPP
import SBayesAPP: option_value, parse_bool, parse_cli_args

function usage()
    println(
        stderr,
        "Usage: julia calculate_multi_chain_genome_summary.jl (--group_path PATH | --chain_paths PATH1,PATH2,...) [--output_path PATH] [--group_label LABEL] [--chain_prefix seed_] [--write_iteration_averages true|false]",
    )
end

function main(args::Vector{String})
    options = parse_cli_args(args; usage=usage, allow_empty=false)

    group_path = option_value(options, "group_path")
    chain_paths = option_value(options, "chain_paths")
    output_path = option_value(options, "output_path")
    group_label = option_value(options, "group_label"; default="Group1")
    chain_prefix = option_value(options, "chain_prefix"; default="seed_")
    write_iteration_averages = parse_bool(option_value(options, "write_iteration_averages"; default="true"))

    if group_path === nothing && chain_paths === nothing
        usage()
        error("Provide either --group_path or --chain_paths")
    end

    result = calculate_multi_chain_genome_summary(
        group_path=group_path,
        chain_paths=chain_paths,
        output_path=output_path,
        chain_prefix=chain_prefix,
        group_label=group_label,
        write_iteration_averages=write_iteration_averages,
    )

    println("Processed chains: ", result.n_chains)
    println("Posterior samples used: ", result.nsamples)
    println("Genome-wide parameter summary written to: ", result.genome_summary_file)
    if !isempty(result.average_file)
        println("Posterior-sample-level genome parameter averages written to: ", result.average_file)
    end
end

try
    main(ARGS)
catch error
    println(stderr, "Error: ", error)
    exit(1)
end
