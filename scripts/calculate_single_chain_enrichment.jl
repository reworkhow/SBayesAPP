using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using SBayesAPP
import SBayesAPP: option_value, parse_bool, parse_cli_args

function usage()
    println(
        stderr,
        "Usage: julia calculate_single_chain_enrichment.jl --data_path PATH --analysis_path PATH --annot_file PATH [--parameter coh2|h21|h22] [--abs_total true|false] [--significance_threshold FLOAT] [--summary_output_file PATH] [--conditional_output_file PATH]",
    )
end

function main(args::Vector{String})
    options = parse_cli_args(args; usage=usage, allow_empty=false)
    data_path = option_value(options, "data_path"; required=true)
    analysis_path = option_value(options, "analysis_path"; required=true)
    annot_file = option_value(options, "annot_file"; required=true)
    parameter = option_value(options, "parameter"; default="coh2")
    abs_total = parse_bool(option_value(options, "abs_total"; default="false"))
    significance_threshold = parse(Float64, option_value(options, "significance_threshold"; default="0.90"))
    summary_output_file = option_value(options, "summary_output_file")
    conditional_output_file = option_value(options, "conditional_output_file")

    result = calculate_single_chain_enrichment(
        data_path,
        analysis_path,
        annot_file;
        parameter=parameter,
        abs_total=abs_total,
        significance_threshold=significance_threshold,
        summary_output_file=summary_output_file,
        conditional_output_file=conditional_output_file,
    )

    println("Posterior enrichment summary written to: ", result.posterior_enrichment_summary_file)
    println("Conditional posterior genetic covariance written to: ", result.conditional_posterior_gcov_file)
end

try
    main(ARGS)
catch error
    println(stderr, "Error: ", error)
    exit(1)
end