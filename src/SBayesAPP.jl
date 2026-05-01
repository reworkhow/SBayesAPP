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
        false,
    )
end

function parse_nonmpi_args(args::Vector{String})
    length(args) == 13 || error(
        "Expected 13 positional arguments for non-MPI mode: data_path analysis_path nIter seed nrank annot_file annot_dict outFreq starting_value_dir secondary_starting_value_dir ST_path thin is_continue",
    )

    return NonMPIConfig(
        _dir_arg(args[1]),
        _dir_arg(args[2]),
        parse(Int, args[3]),
        parse(Int, args[4]),
        parse(Int, args[5]),
        args[6],
        args[7],
        parse(Int, args[8]),
        args[9],
        args[10],
        _dir_arg(args[11]),
        parse(Int, args[12]),
        lowercase(args[13]) == "true",
    )
end

function parse_mpi_args(args::Vector{String})
    length(args) >= 15 || error(
        "Expected at least 15 positional arguments for MPI mode: data_path analysis_path nIter seed nrank annot_file annot_dict outFreq starting_value_dir secondary_starting_value_dir ST_path thin N1 N2 estimate_pi [fixed_hyperparameters] [is_continue] [chr]",
    )

    fixed_hyperparameters = length(args) >= 16 ? lowercase(args[16]) == "true" : false
    is_continue = length(args) >= 17 ? lowercase(args[17]) == "true" : true
    chr = length(args) >= 18 ? args[18] : ""

    return MPIConfig(
        _dir_arg(args[1]),
        _dir_arg(args[2]),
        parse(Int, args[3]),
        parse(Int, args[4]),
        parse(Int, args[5]),
        args[6],
        args[7],
        parse(Int, args[8]),
        args[9],
        args[10],
        _dir_arg(args[11]),
        parse(Int, args[12]),
        parse(Int, args[13]),
        parse(Int, args[14]),
        lowercase(args[15]) == "true",
        fixed_hyperparameters,
        is_continue,
        chr,
    )
end

function build_nonmpi_cmd(config::NonMPIConfig)
    return `$(Base.julia_cmd()) --project=$(repo_root()) $(joinpath(source_root(), "app_nonMPI.jl")) $(config.data_path) $(config.analysis_path) $(string(config.nIter)) $(string(config.seed)) $(string(config.nrank)) $(config.annot_file) $(config.annot_dict) $(string(config.out_freq)) $(config.starting_value_dir) $(config.secondary_starting_value_dir) $(config.st_path) $(string(config.thin)) $(_bool_arg(config.is_continue))`
end

function build_mpi_cmd(config::MPIConfig)
    return `$(Base.julia_cmd()) --project=$(repo_root()) $(joinpath(source_root(), "app_MPI.jl")) $(config.data_path) $(config.analysis_path) $(string(config.nIter)) $(string(config.seed)) $(string(config.nrank)) $(config.annot_file) $(config.annot_dict) $(string(config.out_freq)) $(config.starting_value_dir) $(config.secondary_starting_value_dir) $(config.st_path) $(string(config.thin)) $(string(config.n1)) $(string(config.n2)) $(_bool_arg(config.estimate_pi)) $(_bool_arg(config.fixed_hyperparameters)) $(_bool_arg(config.is_continue)) $(config.chr)`
end

run_nonmpi(config::NonMPIConfig; dry_run::Bool=false) = dry_run ? build_nonmpi_cmd(config) : run_nonmpi_workflow(config)
run_mpi(config::MPIConfig; dry_run::Bool=false) = dry_run ? build_mpi_cmd(config) : run(build_mpi_cmd(config))

end