# Proposal: Refactoring SBayesAPP into a Package-Ready Julia Repository

## Objective

The end goal is not only to optimize the current code, but to turn this repository into a package-ready Julia project that is:

- installable and reproducible,
- organized around reusable modules instead of large scripts,
- easier to test, benchmark, and document,
- prepared for both single-process and parallel execution modes,
- maintainable enough that future optimization does not require editing two giant entrypoint files in parallel.

The current repository already has the key scientific ingredients:

- a non-MPI sampler in `src/app_nonMPI.jl`,
- an MPI sampler in `src/app_MPI.jl`,
- a runnable example under `example/`,
- a code-level description in `code_summary.md`.

What it does not yet have is package structure, dependency management, clean module boundaries, or a refactor plan that separates model logic from execution mode, file I/O, and workflow orchestration.

This proposal focuses on that packaging and refactoring work.

## What package-ready should mean for this repository

For SBayesAPP, package-ready should mean the following.

### 1. Standard Julia project layout

The repository should behave like a normal Julia package, with:

- a `Project.toml` describing dependencies and package metadata,
- an optional `Manifest.toml` for reproducible environments,
- a package entry module under `src/`,
- tests under `test/`,
- example workflows separated from package source code,
- scripts that call package APIs instead of embedding model logic.

### 2. Reusable internal API

The current scientific logic should be callable from Julia functions, not only from top-level scripts driven by `ARGS`.

That means the repository should expose functions such as:

- building a run configuration,
- loading block-wise inputs,
- initializing priors and sampler state,
- running the MCMC sampler,
- writing outputs and restart files,
- launching specific execution modes such as serial or MPI.

### 3. Reproducible environment

A user should be able to clone the repository, activate the Julia environment, instantiate dependencies, and run the example without guessing package versions or relying on a cluster-specific shell environment.

### 4. Clear separation of concerns

The code summary already shows that the current large scripts mix together:

- command-line parsing,
- data loading,
- prior construction,
- sampler state initialization,
- the core marker-update loop,
- posterior summary accumulation,
- checkpoint writing,
- MPI communication.

Those concerns should be separated into modules so that each part can be tested and optimized independently.

## Proposed repository structure

I would restructure the repository toward the following layout.

```text
SBayesAPP/
├── Project.toml
├── Manifest.toml                      # optional to commit, but useful for reproducibility
├── README.md
├── code_summary.md
├── proposal.md
├── docs/                              # optional later phase
├── example/
│   ├── SBayesAPP_input_first10blks/
│   ├── ST_res/
│   ├── annotation_df.txt
│   └── SBayesAPP_res_first10blks/
├── scripts/
│   ├── run_example.jl
│   ├── run_nonmpi.jl
│   └── run_mpi.jl
├── src/
│   ├── SBayesAPP.jl
│   ├── config/
│   │   ├── types.jl
│   │   ├── defaults.jl
│   │   └── cli.jl
│   ├── io/
│   │   ├── inputs.jl
│   │   ├── outputs.jl
│   │   ├── restart.jl
│   │   └── annotations.jl
│   ├── model/
│   │   ├── priors.jl
│   │   ├── state.jl
│   │   ├── statistics.jl
│   │   └── utilities.jl
│   ├── sampler/
│   │   ├── marker_updates.jl
│   │   ├── block_updates.jl
│   │   ├── hyperparameter_updates.jl
│   │   ├── summaries.jl
│   │   └── run_serial.jl
│   ├── parallel/
│   │   ├── mpi.jl
│   │   └── threaded.jl               # later phase
│   └── workflows/
│       ├── nonmpi.jl
│       ├── mpi.jl
│       └── example.jl
├── test/
│   ├── runtests.jl
│   ├── test_config.jl
│   ├── test_priors.jl
│   ├── test_io.jl
│   ├── test_statistics.jl
│   └── test_small_run.jl
└── benchmark/
	├── run_example_benchmark.jl
	└── compare_outputs.jl
```

This structure is intentionally not minimal. It is meant to separate the package core from example data and from user-facing scripts.

## Why this structure fits the current code summary

The existing `code_summary.md` already provides the natural decomposition.

### Configuration layer

The current scripts begin by reading many positional arguments and then converting them into typed values such as:

- paths,
- booleans,
- iteration counts,
- output frequency,
- thinning,
- continuation flags,
- MPI-specific parameters.

This should become a typed configuration layer, for example:

- `RunConfig`
- `SamplerConfig`
- `PathConfig`
- `ParallelConfig`

That change alone would remove a large amount of top-level script noise from both `app_nonMPI.jl` and `app_MPI.jl`.

### Input and output layer

The code summary shows that both scripts spend substantial code on:

- reading annotation tables,
- reading transformed block dictionaries,
- reading single-trait initialization results,
- loading restart files,
- writing posterior summaries,
- writing MCMC samples,
- writing restart outputs.

Those tasks belong in an I/O layer, not inside the sampler loop or entrypoint scripts.

### Model layer

The code summary identifies stable model concepts that should become named types or structured containers:

- `Pi`,
- `A_vec`,
- `R_blk`,
- `betaArray`,
- `deltaArray`,
- `alphaArray`,
- block-wise genetic covariance summaries,
- prior settings such as `Gprior_vec`.

Instead of passing many loose arrays and dictionaries through a giant function, the refactored package should define state structs such as:

- `SamplerState`
- `BlockData`
- `BlockWorkspace`
- `PosteriorAccumulator`
- `PriorState`

These do not need to hide the math. They only need to make ownership and update responsibility explicit.

### Sampler layer

The core MCMC logic described in the code summary is the real heart of the package. That is the code that should eventually be optimized and reused by all execution modes.

From the current summary, the sampler naturally breaks into these units:

- per-marker update logic,
- per-block update logic,
- residual covariance updates,
- `Pi` updates,
- covariance hyperparameter updates,
- posterior summary accumulation,
- checkpoint writing triggers.

Those should each become explicit functions. The current giant loop should become orchestration code that calls them.

### Parallel layer

The MPI version should not duplicate the scientific kernel. It should wrap the shared sampler logic with:

- rank-local data loading,
- reductions to rank 0,
- broadcasts of shared parameters,
- MPI-specific orchestration only.

That means the MPI code should eventually become a thin layer around a common sampler core, not a second implementation of the same algorithm.

## Environment and dependency plan

The repository also needs a clean environment story.

### Immediate environment tasks

I would first create a real Julia `Project.toml` with the currently used dependencies, which appear to include:

- `CSV`
- `DataFrames`
- `Distributions`
- `ProgressMeter`
- `JLD2`
- `MPI` for the MPI workflow only

Standard libraries such as `LinearAlgebra`, `Random`, `Statistics`, `Dates`, and `DelimitedFiles` do not need to be added as external package dependencies, but the package environment should still clearly document what is required.

### Environment policy

I recommend this policy.

- `Project.toml` is committed.
- `Manifest.toml` is committed if exact reproducibility is important for paper or pipeline runs.
- MPI is treated as an optional dependency path in the user documentation, even if the code includes MPI support.
- the README documents one minimal install path for non-MPI users and one MPI path for cluster users.

### Example environment workflow

The intended user workflow should become:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Then either:

```bash
julia scripts/run_example.jl
```

or:

```bash
julia scripts/run_nonmpi.jl --data-path ...
```

with a separate documented MPI workflow for cluster execution.

### Optional future environment improvements

Later, I would also consider:

- a `docs/` environment for documentation builds,
- a `benchmark/Project.toml` if benchmarking grows more complex,
- CI checks for package instantiation and tests.

## How to break down the current giant Julia scripts

The current refactor should be guided directly by the code summary rather than by cosmetic file splitting.

### Current problem

Right now `app_nonMPI.jl` and `app_MPI.jl` each combine:

- argument parsing,
- initialization logic,
- all helper functions,
- file I/O,
- the main MCMC loop,
- summary bookkeeping,
- restart writing,
- execution-mode concerns.

That makes the code hard to:

- test,
- optimize,
- profile,
- reuse,
- package.

### Proposed decomposition of the non-MPI script

I would break `app_nonMPI.jl` into the following logical units.

#### 1. `config/`

Responsibility:

- define typed configs,
- parse CLI arguments,
- validate required paths and options,
- fill defaults.

This replaces most direct `ARGS[...]` access.

#### 2. `io/inputs.jl`

Responsibility:

- load transformed block-wise data,
- load annotation tables,
- load single-trait initialization outputs,
- load block-level sample size inputs,
- normalize file naming conventions.

This moves all file-loading logic out of the main sampler script.

#### 3. `io/restart.jl`

Responsibility:

- load continuation state,
- load last-sample `Pi`, `A_vec`, `betaArray`, `deltaArray`, and `R_blk`,
- write restart outputs in one place.

This removes restart logic from the main MCMC loop and makes continuation behavior testable.

#### 4. `model/priors.jl`

Responsibility:

- build `startPi`,
- build `Gprior_vec`,
- define logic for prior initialization from `ST_path`,
- centralize assumptions about categorical versus continuous annotations.

#### 5. `model/state.jl`

Responsibility:

- define state containers for the sampler,
- initialize arrays for `betaArray`, `deltaArray`, `alphaArray`, `R_blk`, and working buffers,
- expose constructors for fresh runs and continuation runs.

#### 6. `sampler/marker_updates.jl`

Responsibility:

- the per-marker, per-category, per-trait update logic,
- updating `beta`, `delta`, `alpha`, and residual vectors,
- counting mixture states.

This is the most important future optimization target and should be isolated early.

#### 7. `sampler/block_updates.jl`

Responsibility:

- per-block orchestration,
- reconstruction of block-level genetic covariance summaries,
- residual covariance updates.

#### 8. `sampler/hyperparameter_updates.jl`

Responsibility:

- `Pi` updates,
- `A_vec` updates,
- empirical `scale_G_vec` updates,
- related matrix inversions and validation.

#### 9. `sampler/summaries.jl`

Responsibility:

- update posterior means,
- update posterior second moments,
- prepare final summary objects,
- keep all output-accumulation rules in one place.

#### 10. `io/outputs.jl`

Responsibility:

- write final `estA*`, `estB*`, `estG*`, `estPi*`, `estR`, and trace files,
- write MCMC sample files,
- manage output directories and filenames.

#### 11. `workflows/nonmpi.jl`

Responsibility:

- take a typed config,
- call initialization,
- call the shared serial sampler,
- call output writers,
- provide one public function such as `run_nonmpi(config)`.

### Proposed decomposition of the MPI script

The MPI script should be decomposed differently: not by copying the non-MPI code into more files, but by wrapping the same shared kernel with MPI-specific coordination.

I would keep MPI responsibilities limited to:

- MPI initialization and finalization,
- determining `my_rank` and `cluster_size`,
- rank-local input loading,
- reduction of block- and category-level summaries to rank 0,
- broadcast of updated shared parameters back to ranks,
- MPI-specific run orchestration.

This should live under:

- `parallel/mpi.jl`
- `workflows/mpi.jl`

The marker updates, block summaries, hyperparameter updates, and output formatting should not be reimplemented in the MPI file.

## Proposed package API shape

I would aim for a small but clear public API.

At the package level:

- `run_nonmpi(config)`
- `run_mpi(config)`
- `load_example_config()`
- `build_config(; kwargs...)`

Potential internal APIs:

- `load_inputs(config)`
- `initialize_priors(config, inputs)`
- `initialize_state(config, inputs, priors)`
- `run_sampler!(state, inputs, config)`
- `write_outputs(state, accumulators, config)`

This package API would make the code usable from:

- scripts,
- notebooks,
- tests,
- benchmarking harnesses,
- future wrappers.

## Testing strategy for package readiness

A package-ready refactor without tests would still be fragile. I would add tests in layers.

### Unit tests

Test the isolated building blocks:

- argument and config validation,
- prior construction from example ST outputs,
- annotation loading and ordering,
- flatten and unflatten helpers,
- correlation and variance helper functions,
- restart serialization and deserialization.

### Small integration tests

Run a very small serial workflow using the example input or a smaller synthetic fixture to confirm:

- outputs are created,
- shapes are correct,
- continuation mode can reload prior state,
- summary files are internally consistent.

### Regression tests

Use the current example as the initial behavioral baseline.

The goal is not exact chain identity for every stochastic output unless seeds and update order make that possible. The goal is to preserve:

- expected output structure,
- stable interpretation of files,
- reasonable agreement of posterior summaries under controlled settings.

## Script and CLI plan

The repository should still provide easy command-line entrypoints, but these should become thin wrappers around package functions.

### Recommended change

- replace ad hoc shell-first invocation with Julia scripts in `scripts/`,
- keep `script/run.sh` only as a convenience launcher if needed,
- make shell scripts call `julia scripts/run_example.jl` or `julia scripts/run_nonmpi.jl` rather than pointing directly at internal source files.

That keeps the package source clean and makes future changes less brittle.

## Documentation plan

The current repository now has a good start with `README.md` and `code_summary.md`. To support a package-ready repo, documentation should be split by purpose.

### README

Should answer:

- what SBayesAPP does,
- how to install the environment,
- how to run the example,
- where outputs go,
- where to find detailed code and method notes.

### `code_summary.md`

Should remain the internal developer-oriented description of the current algorithm and current execution flow.

### Future docs

I would eventually add:

- input format documentation,
- output file documentation,
- continuation workflow documentation,
- package API usage examples,
- parallel execution notes.

## Phased implementation plan

I would implement the full refactor in phases.

### Phase 1: Make the repository a real Julia project

Deliverables:

- `Project.toml`
- package module entry file `src/SBayesAPP.jl`
- initial README environment instructions
- scripts moved to a stable `scripts/` directory

This phase makes the repository installable before the deeper refactor starts.

### Phase 2: Extract non-MPI code into reusable modules

Deliverables:

- typed configuration objects,
- input/output modules,
- prior and state modules,
- serial workflow function,
- non-MPI script reduced to a thin wrapper.

This phase creates the package core.

### Phase 3: Refactor MPI around the shared kernel

Deliverables:

- MPI-specific orchestration layer,
- shared sampler logic reused from the serial path,
- clearer separation between communication and model updates.

This phase prevents long-term divergence between serial and parallel code.

### Phase 4: Add tests, benchmarks, and regression harnesses

Deliverables:

- `test/` suite,
- example benchmark script,
- baseline output comparison tooling.

This phase makes later optimization safe.

### Phase 5: Optimize the refactored sampler core

Only after the package structure and tests exist should the main performance work begin.

Targets:

- inner marker-update loop,
- allocation reduction,
- summary accumulation,
- restart I/O efficiency,
- potential threaded execution path.

### Phase 6: Add a user-friendly shared-memory parallel path

At this phase, I would consider a threaded serial-core variant as the first easier-to-use alternative to MPI.

I would not start with this before the package refactor, because parallelizing unstable large scripts would make the code harder to clean up later.

## Recommendation

The right way to package this repository is:

1. convert it into a real Julia project,
2. factor the current scripts into modules based on the boundaries already identified in `code_summary.md`,
3. make non-MPI the clean reference workflow,
4. rebuild MPI as a thin distributed layer around the shared kernel,
5. add tests and benchmarks before serious optimization,
6. then optimize and parallelize the refactored core.

This path aligns with the actual structure of the current code and gives a realistic route from a working research repository to a package-ready scientific software project.