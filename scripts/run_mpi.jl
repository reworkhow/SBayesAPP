using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using SBayesAPP

run_mpi(parse_mpi_args(ARGS))