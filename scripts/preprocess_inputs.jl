using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using SBayesAPP
import SBayesAPP: option_value, parse_cli_args

function usage()
    println(
        stderr,
        "Usage: julia --project=. scripts/preprocess_inputs.jl " *
        "--LD_info_path PATH --out PATH --trait1_file FILE --trait2_file FILE --LDinfo_file FILE --annot_file FILE [options]\n\n" *
        "Options:\n" *
        "  --annot_dict_name NAME       Annotation dict filename stem [default: anno_matrix_dict]\n" *
        "  --readable_files_dir NAME    LD readable files subdirectory [default: readableFiles]\n" *
        "  --nblock INT                 Process only the first n blocks from LD info [default: all]\n"
    )
end

function main(args::Vector{String})
    options = parse_cli_args(args; usage=usage, allow_empty=false)

    result = build_nonmpi_input_dicts(
        option_value(options, "out"; required=true),
        option_value(options, "LD_info_path"; required=true),
        option_value(options, "trait1_file"; required=true),
        option_value(options, "trait2_file"; required=true),
        option_value(options, "LDinfo_file"; required=true),
        option_value(options, "annot_file"; required=true);
        annot_dict_name=option_value(options, "annot_dict_name"; default="anno_matrix_dict"),
        readable_files_dir=option_value(options, "readable_files_dir"; default="readableFiles"),
        nblocks=let value = option_value(options, "nblock")
            isnothing(value) ? nothing : parse(Int, value)
        end,
    )

    println("Preprocess completed.")
    println("output_path: " * result.output_path)
    println("annotation_file: " * result.annotation_file)
    println("annot_dict: " * result.annot_dict)
    println("annotation_dict_file: " * result.annotation_dict_file)
end

try
    main(ARGS)
catch error
    println(stderr, "Error: ", error)
    exit(1)
end