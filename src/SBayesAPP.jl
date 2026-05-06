module SBayesAPP

include("config/types.jl")
include("io/annotations.jl")
include("io/inputs.jl")
include("io/restart.jl")
include("model/continuation_state.jl")
include("model/mcmc_setup.jl")
include("model/priors.jl")
include("model/r_blk_state.jl")
include("model/state.jl")
include("model/block_setup.jl")
include("model/utilities.jl")
include("io/outputs.jl")
include("workflows/nonmpi.jl")

using .ConfigTypes: MPIConfig, NonMPIConfig 

export MPIConfig,
       NonMPIConfig,
       build_mpi_cmd,
       build_nonmpi_cmd,
    build_gprior_vec,
    build_start_pi,
       example_nonmpi_config,
    flatten,
    flatten_matrices,
    compute_correlation,
    is_positive_definite,
    load_annotation_metadata,
     load_nonmpi_block_data,
       parse_mpi_args,
       parse_nonmpi_args,
    read_to_dict,
    read_to_dict_posterior_mean,
       repo_root,
       run_mpi,
       run_nonmpi,
    source_root,
    unflatten,
    unflatten_matrices

repo_root() = normpath(joinpath(@__DIR__, ".."))
source_root() = joinpath(repo_root(), "src")

_bool_arg(value::Bool) = value ? "true" : "false"
_dir_arg(path::AbstractString) = endswith(path, "/") ? String(path) : string(path, "/")
_normalize_cli_key(key::AbstractString) = lowercase(replace(String(key), "-" => "_"))

function _parse_bool_arg(value::AbstractString)
    lower_value = lowercase(value)
    lower_value == "true" && return true
    lower_value == "false" && return false
    error("Expected boolean argument 'true' or 'false', got: $value")
end

function _parse_int_arg(value::AbstractString, option_name::AbstractString)
    try
        return parse(Int, value)
    catch error
        throw(ArgumentError("Expected integer value for --$option_name, got: $value"))
    end
end

function _parse_named_cli_args(args::Vector{String})
    parsed = Dict{String,String}()
    index = 1
    while index <= length(args)
        token = args[index]
        startswith(token, "--") || error("Expected named option starting with --, got: $token")

        option = token[3:end]
        value = nothing
        if occursin('=', option)
            option, value = split(option, '='; limit=2)
        else
            index == length(args) && error("Missing value for option --$option")
            value = args[index + 1]
            index += 1
        end

        normalized_option = _normalize_cli_key(option)
        haskey(parsed, normalized_option) && error("Duplicate option --$option")
        parsed[normalized_option] = value
        index += 1
    end
    return parsed
end

function _lookup_named_arg(parsed::Dict{String,String}, option_names::Vector{String}; required::Bool=false, default=nothing)
    for option_name in option_names
        normalized_option = _normalize_cli_key(option_name)
        if haskey(parsed, normalized_option)
            return parsed[normalized_option]
        end
    end

    required && error("Missing required option --$(first(option_names))")
    return default
end

_has_named_cli_args(args::Vector{String}) = any(startswith(arg, "--") for arg in args)

_default_nonmpi_cli_settings() = (
    n_con=0,
    estimate_vare=true,
    estimate_vara=true,
    estimate_pi=true,
    estimate_Gscale=true,
    estGscale_iter=500,
    report_pleiotropic_qtl_effect_matrix=false,
    output_mcmc_delta=false,
)

function example_nonmpi_config(; root::AbstractString=repo_root())
    return NonMPIConfig(
        _dir_arg(joinpath(root, "example", "SBayesAPP_input_first10blks")),
        _dir_arg(joinpath(root, "example", "SBayesAPP_res_first10blks")),
        1000,
        42,
        1,
        "annotation_df.txt",
        "anno_matrix_dict",
        100,
        "XXX",
        "XXX",
        _dir_arg(joinpath(root, "example", "ST_res")),
        50,
        300000,
        300000,
        false,
        n_con=0,
        estimate_vare=true,
        estimate_vara=true,
        estimate_pi=true,
        estimate_Gscale=true,
        estGscale_iter=500,
        report_pleiotropic_qtl_effect_matrix=true,
        output_mcmc_delta=false,
    )
end

function parse_nonmpi_args(args::Vector{String})
    _has_named_cli_args(args) || error(
        "Expected named options for non-MPI mode, for example --data_path ... --analysis_path ...",
    )

    parsed = _parse_named_cli_args(args)
    defaults = _default_nonmpi_cli_settings()

    return NonMPIConfig(
        _dir_arg(_lookup_named_arg(parsed, ["data_path"]; required=true)),
        _dir_arg(_lookup_named_arg(parsed, ["analysis_path"]; required=true)),
        _parse_int_arg(_lookup_named_arg(parsed, ["n_iter", "niter"]; required=true), "n_iter"),
        _parse_int_arg(_lookup_named_arg(parsed, ["seed"]; required=true), "seed"),
        _parse_int_arg(_lookup_named_arg(parsed, ["nrank"]; required=true), "nrank"),
        _lookup_named_arg(parsed, ["annot_file"]; required=true),
        _lookup_named_arg(parsed, ["annot_dict"]; required=true),
        _parse_int_arg(_lookup_named_arg(parsed, ["out_freq", "outfreq"]; required=true), "out_freq"),
        _lookup_named_arg(parsed, ["starting_value_dir"]; required=true),
        _lookup_named_arg(parsed, ["gscale_value_dir", "secondary_starting_value_dir"]; required=true),
        _dir_arg(_lookup_named_arg(parsed, ["st_path"]; required=true)),
        _parse_int_arg(_lookup_named_arg(parsed, ["thin"]; required=true), "thin"),
        _parse_int_arg(_lookup_named_arg(parsed, ["n1"]; required=true), "n1"),
        _parse_int_arg(_lookup_named_arg(parsed, ["n2"]; required=true), "n2"),
        _parse_bool_arg(_lookup_named_arg(parsed, ["is_continue"]; required=true)),
        n_con=_parse_int_arg(string(_lookup_named_arg(parsed, ["n_con", "ncon"]; default=string(defaults.n_con))), "n_con"),
        estimate_vare=_parse_bool_arg(string(_lookup_named_arg(parsed, ["estimate_vare"]; default=string(defaults.estimate_vare)))),
        estimate_vara=_parse_bool_arg(string(_lookup_named_arg(parsed, ["estimate_vara"]; default=string(defaults.estimate_vara)))),
        estimate_pi=_parse_bool_arg(string(_lookup_named_arg(parsed, ["estimate_pi"]; default=string(defaults.estimate_pi)))),
        estimate_Gscale=_parse_bool_arg(string(_lookup_named_arg(parsed, ["estimate_gscale"]; default=string(defaults.estimate_Gscale)))),
        estGscale_iter=_parse_int_arg(string(_lookup_named_arg(parsed, ["estgscale_iter"]; default=string(defaults.estGscale_iter))), "estGscale_iter"),
        report_pleiotropic_qtl_effect_matrix=_parse_bool_arg(string(_lookup_named_arg(parsed, ["report_pleiotropic_qtl_effect_matrix"]; default=string(defaults.report_pleiotropic_qtl_effect_matrix)))),
        output_mcmc_delta=_parse_bool_arg(string(_lookup_named_arg(parsed, ["output_mcmc_delta"]; default=string(defaults.output_mcmc_delta)))),
    )
end

function parse_mpi_args(args::Vector{String})
    _has_named_cli_args(args) || error(
        "Expected named options for MPI mode, for example --data_path ... --analysis_path ...",
    )

    parsed = _parse_named_cli_args(args)

    return MPIConfig(
        _dir_arg(_lookup_named_arg(parsed, ["data_path"]; required=true)),
        _dir_arg(_lookup_named_arg(parsed, ["analysis_path"]; required=true)),
        _parse_int_arg(_lookup_named_arg(parsed, ["n_iter", "niter"]; required=true), "n_iter"),
        _parse_int_arg(_lookup_named_arg(parsed, ["seed"]; required=true), "seed"),
        _parse_int_arg(_lookup_named_arg(parsed, ["nrank"]; required=true), "nrank"),
        _lookup_named_arg(parsed, ["annot_file"]; required=true),
        _lookup_named_arg(parsed, ["annot_dict"]; required=true),
        _parse_int_arg(_lookup_named_arg(parsed, ["out_freq", "outfreq"]; required=true), "out_freq"),
        _lookup_named_arg(parsed, ["starting_value_dir"]; required=true),
        _lookup_named_arg(parsed, ["secondary_starting_value_dir"]; required=true),
        _dir_arg(_lookup_named_arg(parsed, ["st_path"]; required=true)),
        _parse_int_arg(_lookup_named_arg(parsed, ["thin"]; required=true), "thin"),
        _parse_int_arg(_lookup_named_arg(parsed, ["n1"]; required=true), "n1"),
        _parse_int_arg(_lookup_named_arg(parsed, ["n2"]; required=true), "n2"),
        _parse_int_arg(string(_lookup_named_arg(parsed, ["n_con", "ncon"]; default="0")), "n_con"),
        _parse_bool_arg(_lookup_named_arg(parsed, ["estimate_pi"]; required=true)),
        _parse_bool_arg(string(_lookup_named_arg(parsed, ["fixed_hyperparameters"]; default="false"))),
        _parse_bool_arg(string(_lookup_named_arg(parsed, ["is_continue"]; default="true"))),
        string(_lookup_named_arg(parsed, ["chr"]; default="")),
    )
end

function build_nonmpi_cmd(config::NonMPIConfig)
    return `$(Base.julia_cmd()) --project=$(repo_root()) $(joinpath(source_root(), "app_nonMPI.jl")) --data_path $(config.data_path) --analysis_path $(config.analysis_path) --n_iter $(string(config.nIter)) --seed $(string(config.seed)) --nrank $(string(config.nrank)) --annot_file $(config.annot_file) --annot_dict $(config.annot_dict) --out_freq $(string(config.out_freq)) --starting_value_dir $(config.starting_value_dir) --gscale_value_dir $(config.gscale_value_dir) --st_path $(config.st_path) --thin $(string(config.thin)) --n1 $(string(config.n1)) --n2 $(string(config.n2)) --n_con $(string(config.n_con)) --is_continue $(_bool_arg(config.is_continue)) --estimate_vare $(_bool_arg(config.estimate_vare)) --estimate_vara $(_bool_arg(config.estimate_vara)) --estimate_pi $(_bool_arg(config.estimate_pi)) --estimate_gscale $(_bool_arg(config.estimate_Gscale)) --estgscale_iter $(string(config.estGscale_iter)) --report_pleiotropic_qtl_effect_matrix $(_bool_arg(config.report_pleiotropic_qtl_effect_matrix)) --output_mcmc_delta $(_bool_arg(config.output_mcmc_delta))`
end

function build_mpi_cmd(config::MPIConfig)
    return `$(Base.julia_cmd()) --project=$(repo_root()) $(joinpath(source_root(), "app_MPI.jl")) --data_path $(config.data_path) --analysis_path $(config.analysis_path) --n_iter $(string(config.nIter)) --seed $(string(config.seed)) --nrank $(string(config.nrank)) --annot_file $(config.annot_file) --annot_dict $(config.annot_dict) --out_freq $(string(config.out_freq)) --starting_value_dir $(config.starting_value_dir) --secondary_starting_value_dir $(config.secondary_starting_value_dir) --st_path $(config.st_path) --thin $(string(config.thin)) --n1 $(string(config.n1)) --n2 $(string(config.n2)) --n_con $(string(config.n_con)) --estimate_pi $(_bool_arg(config.estimate_pi)) --fixed_hyperparameters $(_bool_arg(config.fixed_hyperparameters)) --is_continue $(_bool_arg(config.is_continue)) --chr $(config.chr)`
end

run_nonmpi(config::NonMPIConfig; dry_run::Bool=false) = dry_run ? build_nonmpi_cmd(config) : run_nonmpi_workflow(config)
run_mpi(config::MPIConfig; dry_run::Bool=false) = dry_run ? build_mpi_cmd(config) : run(build_mpi_cmd(config))

end