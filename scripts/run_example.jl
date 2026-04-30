using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using SBayesAPP

run_nonmpi(example_nonmpi_config())