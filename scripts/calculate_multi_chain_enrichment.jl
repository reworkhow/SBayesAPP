using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using SBayesAPP
import SBayesAPP: option_value, parse_bool, parse_cli_args

function usage()
    println(
        stderr,
        "Usage: julia calculate_multi_chain_enrichment.jl --data_path PATH --annot_file PATH (--group_path PATH | --chain_paths PATH1,PATH2,...) [--output_path PATH] [--group_label LABEL] [--chain_prefix seed_] [--enrichment_parameter coh2|h21|h22] [--abs_total true|false] [--significance_threshold FLOAT] [--write_iteration_averages true|false]",
    )
end


function main(args::Vector{String})
    options = parse_cli_args(args; usage=usage, allow_empty=false)

    data_path = option_value(options, "data_path"; required=true)
    annot_file = option_value(options, "annot_file"; required=true)
    group_path = option_value(options, "group_path")
    chain_paths = option_value(options, "chain_paths")
    output_path = option_value(options, "output_path")
    group_label = option_value(options, "group_label"; default="Group1")
    chain_prefix = option_value(options, "chain_prefix"; default="seed_")
    enrichment_parameter = option_value(options, "enrichment_parameter"; default="coh2")
    abs_total = parse_bool(option_value(options, "abs_total"; default="false"))
    significance_threshold = parse(Float64, option_value(options, "significance_threshold"; default="0.90"))
    write_iteration_averages = parse_bool(option_value(options, "write_iteration_averages"; default="true"))

    if group_path === nothing && chain_paths === nothing
        usage()
        error("Provide either --group_path or --chain_paths") 
    end

    result = calculate_multi_chain_enrichment(
        data_path,
        annot_file;
        group_path=group_path,
        chain_paths=chain_paths,
        output_path=output_path,
        chain_prefix=chain_prefix,
        group_label=group_label,
        enrichment_parameter=enrichment_parameter,
        abs_total=abs_total,
        significance_threshold=significance_threshold,
        write_iteration_averages=write_iteration_averages,
    )

    println("Processed chains: ", result.n_chains)
    println("Posterior samples used: ", result.nsamples)
    println("Enrichment summary written to: ", result.enrichment_summary_file)
    println("Conditional posterior genetic covariance written to: ", result.conditional_summary_file)
    println("SNP effect correlation summary written to: ", result.snpcor_summary_file)
    println("Pi/polygenicity summary written to: ", result.pi_summary_file)
    if !isempty(result.average_files)
        println("Posterior-sample-level group averages written to: ", join(result.average_files, ", "))
    end
end

try
    main(ARGS)
catch error
    println(stderr, "Error: ", error)
    exit(1)
end
