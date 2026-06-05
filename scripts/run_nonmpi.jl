using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using SBayesAPP

run_nonmpi(parse_nonmpi_args(ARGS))