# Proposal: Refactoring SBayesAPP into a Package-Ready, Threaded Julia Repository

## Objective

The end goal is not only to optimize the current code, but to turn this repository into a package-ready Julia project that is:

- installable and reproducible,
- organized around reusable modules instead of large scripts,
- easier to test, benchmark, and document,
- prepared for both serial and shared-memory parallel execution,
- maintainable enough that future optimization does not require carrying a second distributed implementation in parallel.

The repository already has the key scientific ingredients:

- a package-backed non-MPI workflow centered on `src/workflows/nonmpi.jl`,
- command-line wrappers under `scripts/`,
- a runnable example under `example/`,
- tests under `test/`,
- a code-level description in `code_summary.md`.

What it still needs is a refactor plan that treats the current non-MPI workflow as the single reference implementation and adds Julia thread-based parallelism around that shared kernel instead of maintaining an MPI-specific code path.

This proposal focuses on that packaging and threading work.

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

The scientific logic should be callable from Julia functions, not only from top-level scripts driven by `ARGS`.

That means the repository should expose functions such as:

- building a run configuration,
- loading block-wise inputs,
- initializing priors and sampler state,
- running the MCMC sampler,
- writing outputs and restart files,
- launching either a serial or threaded execution policy over the same sampler core.

### 3. Reproducible environment

A user should be able to clone the repository, activate the Julia environment, instantiate dependencies, and run the example without guessing package versions or relying on a cluster-specific shell environment.

### 4. Clear separation of concerns

The current workflow still mixes together:

- command-line parsing,
- data loading,
- prior construction,
- sampler state initialization,
- the core marker-update loop,
- posterior summary accumulation,
- checkpoint writing,
- and future parallel execution concerns.

Those concerns should be separated into modules so that each part can be tested and optimized independently.

### 5. Thread-safe parallel design

Because the intended parallel strategy is Julia multi-threading rather than MPI, package-ready for this repository must also mean:

- no hidden data races in the sampler core,
- no dependence on mutating shared `Dict` or `Vector` state from multiple tasks without a clear ownership rule,
- explicit local reductions for quantities that are accumulated across blocks,
- a serial fallback path that remains the reference for correctness.

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
│   └── run_nonmpi.jl
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
│   │   ├── utilities.jl
│   │   └── block_data.jl
│   ├── sampler/
│   │   ├── marker_updates.jl
│   │   ├── block_updates.jl
│   │   ├── hyperparameter_updates.jl
│   │   ├── summaries.jl
│   │   ├── run_serial.jl
│   │   └── reductions.jl
│   ├── parallel/
│   │   └── threaded.jl
│   └── workflows/
│       ├── nonmpi.jl
│       └── example.jl
├── test/
│   ├── runtests.jl
│   ├── test_config.jl
│   ├── test_priors.jl
│   ├── test_io.jl
│   ├── test_statistics.jl
│   ├── test_small_run.jl
│   └── test_threaded_run.jl
└── benchmark/
	├── run_example_benchmark.jl
	├── compare_outputs.jl
	└── benchmark_threaded_scaling.jl
```

This structure is intentionally not minimal. It is meant to separate the package core from example data and from user-facing scripts while reserving one small area for threading-specific orchestration.

## Why this structure fits the current code

The current code already suggests the right decomposition.

### Configuration layer

The current workflow reads many named arguments and converts them into typed values such as:

- paths,
- booleans,
- iteration counts,
- output frequency,
- thinning,
- continuation flags,
- output toggles,
- and future execution policy settings.

This should become a typed configuration layer, for example:

- `RunConfig`
- `SamplerConfig`
- `PathConfig`
- `ParallelConfig`

For the threaded path, `ParallelConfig` should describe behavior such as:

- `mode = :serial | :threads`,
- block chunk size,
- scheduling policy,
- whether block updates are threaded,
- whether certain summaries remain serial.

The number of Julia threads itself should not be treated as a normal runtime argument, because Julia thread count is configured before startup via `--threads` or `JULIA_NUM_THREADS`.

### Input and output layer

The current workflow spends substantial code on:

- reading annotation tables,
- reading transformed block dictionaries,
- reading single-trait initialization results,
- loading restart files,
- writing posterior summaries,
- writing MCMC samples,
- writing restart outputs.

Those tasks belong in an I/O layer, not inside the sampler loop or entrypoint scripts.

### Model layer

The current code has stable model concepts that should become named types or structured containers:

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
- `ThreadLocalAccumulator`

These do not need to hide the math. They only need to make ownership and update responsibility explicit.

### Sampler layer

The core MCMC logic is the real heart of the package. That is the code that should be optimized and reused by both serial and threaded execution.

From the current workflow, the sampler naturally breaks into these units:

- per-marker update logic,
- per-block update logic,
- residual covariance updates,
- `Pi` updates,
- covariance hyperparameter updates,
- posterior summary accumulation,
- checkpoint writing triggers.

Those should each become explicit functions. The current large loop should become orchestration code that calls them.

### Parallel layer

The threaded version should not become a second implementation of the model. It should wrap the shared sampler logic with:

- block partitioning,
- per-task local accumulators,
- explicit reduction of block-level summaries,
- thread-safe orchestration only.

That means the threaded code should be a thin layer around a common sampler core, not a fork of the model implementation.

## Thread-based execution plan

This is the key design shift relative to the older MPI-oriented plan.

### Guiding principle

The correct reference behavior should remain the serial non-MPI workflow. Threading should be added by parallelizing independent block work inside each iteration, while keeping global hyperparameter updates and final reductions explicit and testable.

### Recommended parallel unit: LD blocks

The best first threading target is the block loop, not the inner marker loop.

That is the right first cut because:

- LD blocks are already the natural decomposition of the data,
- each block owns its own transformed design and response state,
- block-local updates can be run with shared global hyperparameters but independent local residual and effect updates,
- the resulting quantities can be reduced after the threaded region.

In other words, the intended iteration shape should become:

1. freeze global parameters for the current iteration,
2. process blocks in parallel,
3. return local block summaries from each task,
4. reduce those summaries,
5. update global hyperparameters serially,
6. write outputs and checkpoints as before.

### Data layout changes needed before threading

Before threading the block loop, the package should stop relying on mutable `Dict` objects inside hot parallel regions.

In particular, block-aligned vectors or structs should replace repeated keyed access such as:

- transformed `Y` by block,
- annotation masks by block,
- SNP indices by block,
- per-block sample sizes,
- per-block working buffers.

This is important for two reasons:

- it reduces lookup overhead in the hot path,
- it avoids unsafe concurrent mutation of base collections like `Dict`.

The package should construct a `Vector{BlockData}` or similar once during setup, where each entry contains the full state for one block.

### Local accumulation strategy

Julia threading does not provide an automatic reduction argument for `Threads.@threads`, and shared mutation is the main correctness risk. Therefore, the threaded design should use per-task local accumulators rather than shared counters.

Examples of quantities that should be accumulated locally and then reduced include:

- `nLoci_array_vec`,
- `SSE_vec`,
- any block-level summaries that currently aggregate into global counters,
- temporary genetic covariance summaries when they are not stored in explicitly block-owned slots.

This suggests a design where each task returns a small reduction object, for example:

- `ThreadLocalAccumulator`
- `ThreadLocalHyperStats`
- `ThreadLocalOutputState`

These objects can then be merged serially after `fetch`.

### Safe shared writes

Some writes can remain in place during threaded execution, but only if ownership is explicit.

Safe examples include writing to disjoint block-owned or SNP-slice-owned storage such as:

- `R_blk[b]`,
- block-indexed variance arrays like `varg_blk_cat[b]`,
- disjoint SNP slices of `betaArray`, `alphaArray`, and `deltaArray` corresponding to non-overlapping `SNPIndexb`.

Unsafe examples include:

- mutating a shared `Dict` from multiple tasks,
- `push!` into shared arrays,
- using one global temporary buffer across tasks,
- relying on `threadid()` to choose a buffer for a yielding task.

### Scheduling recommendation

The first implementation should prefer task-based chunking with `Threads.@spawn` over a naive `Threads.@threads` rewrite of the full loop.

That is the safer route because:

- block sizes are heterogeneous,
- chunked tasks make local reductions straightforward,
- task-returned results are easier to reason about than shared mutable state,
- Julia task migration makes `threadid()`-indexed buffers fragile.

The intended pattern is:

1. partition blocks into chunks,
2. `Threads.@spawn` one task per chunk,
3. allocate local accumulators inside the task,
4. process the assigned blocks,
5. return local summaries,
6. reduce on the caller.

Only once that path is correct should finer-grained scheduling experiments be considered.

### Threading constraints from Julia's model

The official Julia threading model implies several design constraints that should be made explicit in the package plan.

- Julia thread count is configured before startup using `--threads` or `JULIA_NUM_THREADS`; it is not a normal package parameter.
- `Threads.@threads` does not provide a built-in reduction clause, so reductions must be implemented explicitly.
- Task migration means code should not rely on `threadid()` being stable across yields.
- Base collections such as `Dict` require manual locking if they are modified from multiple threads.
- Side-effect-heavy code must be audited carefully before being moved into threaded regions.
- Locks or atomics should be a last resort in the hot path; per-task local state plus serial reduction is preferable for most sampler summaries.

These constraints argue for a reduction-oriented block scheduler rather than shared-state threading.

## Environment and dependency plan

The repository also needs a clean environment story.

### Immediate environment tasks

The Julia environment should include the packages already needed by the current non-MPI path, such as:

- `CSV`
- `DataFrames`
- `Distributions`
- `ProgressMeter`
- `JLD2`

Standard libraries such as `LinearAlgebra`, `Random`, `Statistics`, `Dates`, `DelimitedFiles`, and `Base.Threads` do not need to be added as external dependencies.

Unlike MPI, Julia threads do not require a separate communication package for the initial implementation.

### Environment policy

I recommend this policy.

- `Project.toml` is committed.
- `Manifest.toml` is committed if exact reproducibility is important for paper or pipeline runs.
- Threaded execution is documented as the preferred acceleration path.
- README examples show both a normal serial launch and a threaded launch using Julia's startup flags.

### Example environment workflow

The intended user workflow should become:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Then either:

```bash
bash script/run.sh
```

or:

```bash
julia --project=. --threads auto scripts/run_nonmpi.jl --data_path ...
```

### Optional future environment improvements

Later, I would also consider:

- a `docs/` environment for documentation builds,
- a `benchmark/Project.toml` if benchmarking grows more complex,
- CI checks for package instantiation and tests,
- documented guidance for `--gcthreads` once the threaded sampler starts allocating heavily enough that GC tuning matters.

## How to break down the current workflow code

The current refactor should be guided directly by the model structure rather than by cosmetic file splitting.

### Current problem

Right now the non-MPI path still combines:

- argument parsing,
- initialization logic,
- helper functions,
- file I/O,
- the main MCMC loop,
- summary bookkeeping,
- restart writing,
- and future parallel-execution concerns.

That makes the code hard to:

- test,
- optimize,
- profile,
- reuse,
- package,
- and thread safely.

### Proposed decomposition of the non-MPI workflow

I would break the current path into the following logical units.

#### 1. `config/`

Responsibility:

- define typed configs,
- parse CLI arguments,
- validate required paths and options,
- fill defaults.

This replaces direct `ARGS[...]` access and keeps future execution-policy settings centralized.

#### 2. `io/inputs.jl`

Responsibility:

- load transformed block-wise data,
- load annotation tables,
- load single-trait initialization outputs,
- load block-level sample size inputs,
- normalize file naming conventions.

This moves all file-loading logic out of the main sampler workflow.

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

#### 6. `model/block_data.jl`

Responsibility:

- transform block dictionaries into a block-aligned vector representation,
- make block ownership explicit,
- prepare thread-safe read-mostly inputs before entering the sampler.

This is a key prerequisite for threading.

#### 7. `sampler/marker_updates.jl`

Responsibility:

- the per-marker, per-category, per-trait update logic,
- updating `beta`, `delta`, `alpha`, and block-local residual vectors,
- counting mixture states for the current local accumulator.

This is the most important future optimization target and should be isolated early.

#### 8. `sampler/block_updates.jl`

Responsibility:

- per-block orchestration,
- reconstruction of block-level genetic covariance summaries,
- residual covariance updates,
- returning block-local reduction statistics.

#### 9. `sampler/hyperparameter_updates.jl`

Responsibility:

- `Pi` updates,
- `A_vec` updates,
- empirical `scale_G_vec` updates,
- related matrix inversions and validation.

These remain serial until there is a compelling reason to parallelize them.

#### 10. `sampler/summaries.jl`

Responsibility:

- update posterior means,
- update posterior second moments,
- prepare final summary objects,
- keep all output-accumulation rules in one place.

#### 11. `sampler/reductions.jl`

Responsibility:

- define merge rules for thread-local accumulators,
- keep serial and threaded reductions consistent,
- make correctness checks and regression tests easier.

#### 12. `parallel/threaded.jl`

Responsibility:

- partition blocks into chunks,
- spawn tasks,
- allocate per-task workspaces and local accumulators,
- collect and reduce task results,
- provide one thread-aware execution policy without changing the model code.

#### 13. `workflows/nonmpi.jl`

Responsibility:

- take a typed config,
- call initialization,
- choose a serial or threaded block execution policy,
- call output writers,
- provide one public function such as `run_nonmpi(config)`.

## Proposed package API shape

I would aim for a small but clear public API.

At the package level:

- `run_nonmpi(config)`
- `example_nonmpi_config()`
- `build_config(; kwargs...)`

Optionally, after the threaded path exists:

- `run_nonmpi(config; execution=:serial)`
- `run_nonmpi(config; execution=:threads)`

Potential internal APIs:

- `load_inputs(config)`
- `initialize_priors(config, inputs)`
- `initialize_state(config, inputs, priors)`
- `run_iteration_serial!(state, inputs, config)`
- `run_iteration_threaded!(state, inputs, config)`
- `reduce_thread_locals!(...)`
- `write_outputs(state, accumulators, config)`

This package API would make the code usable from:

- scripts,
- notebooks,
- tests,
- benchmarking harnesses,
- and future wrappers.

## Testing strategy for package readiness and threading

A package-ready refactor without tests would still be fragile. Threading raises that bar further.

### Unit tests

Test the isolated building blocks:

- argument and config validation,
- prior construction from example ST outputs,
- annotation loading and ordering,
- flatten and unflatten helpers,
- correlation and variance helper functions,
- restart serialization and deserialization,
- local reduction merge logic.

### Small integration tests

Run a very small workflow using the example input or a smaller synthetic fixture to confirm:

- outputs are created,
- shapes are correct,
- continuation mode can reload prior state,
- summary files are internally consistent.

### Serial versus threaded regression tests

This is the most important new test category.

The threaded path should be checked against the serial path on a fixed small run for:

- matching output structure,
- agreement of posterior summaries within a documented tolerance,
- consistent restart outputs,
- correct behavior for `Threads.nthreads() == 1`, `2`, and a larger count.

The goal is not necessarily bitwise equality once scheduling changes, but it should be possible to verify that:

- the same algorithm is being run,
- reductions are correct,
- no thread-specific file or state corruption appears,
- outputs remain statistically compatible.

### Benchmark tests

Benchmarking should measure:

- total runtime,
- allocation count,
- scaling with thread count,
- whether thread overhead dominates on small examples,
- and whether block-size imbalance hurts scaling.

## Script and CLI plan

The repository should still provide easy command-line entrypoints, but these should remain thin wrappers around package functions.

### Recommended change

- keep `scripts/run_nonmpi.jl` as the main entrypoint,
- keep `script/run.sh` only as a convenience launcher if desired,
- document threaded runs by changing how Julia is launched, not by reviving a separate MPI script.

That means the recommended threaded invocation becomes:

```bash
julia --project=. --threads auto scripts/run_nonmpi.jl --data_path ...
```

If a convenience wrapper is desired later, it should only validate that Julia was started with more than one thread. It should not attempt to create threads itself, because thread count is fixed at process startup.

## Documentation plan

The current repository already has a good start with `README.md` and `code_summary.md`. To support a package-ready threaded repo, documentation should be split by purpose.

### README

Should answer:

- what SBayesAPP does,
- how to install the environment,
- how to run the example,
- how to run with `--threads auto`,
- where outputs go,
- where to find detailed code and method notes.

### `code_summary.md`

Should remain the internal developer-oriented description of the current algorithm and current execution flow.

It should eventually include a short note describing:

- which iteration stages are block-parallel,
- which accumulators are reduced after the threaded region,
- which steps remain intentionally serial.

### Future docs

I would eventually add:

- input format documentation,
- output file documentation,
- continuation workflow documentation,
- package API usage examples,
- threaded execution notes,
- performance tuning notes for thread count and chunk size.

## Phased implementation plan

I would implement the refactor in the following phases.

### Phase 1: Consolidate the repository around the current package-backed non-MPI path

Deliverables:

- one clear package module entry file `src/SBayesAPP.jl`,
- one primary CLI wrapper in `scripts/run_nonmpi.jl`,
- README wording aligned to the non-MPI path,
- removal of any remaining MPI-first design assumptions from planning and docs.

This phase makes the non-MPI workflow the single reference implementation.

### Phase 2: Extract the serial sampler into reusable modules

Deliverables:

- typed configuration objects,
- input/output modules,
- prior and state modules,
- a serial workflow function,
- block-aligned input/state containers,
- the current non-MPI script reduced to a thin wrapper.

This phase creates the package core.

### Phase 3: Make the serial core thread-ready

Deliverables:

- no mutable `Dict` dependence in hot loop ownership,
- explicit per-block state containers,
- explicit local reduction objects,
- serial helper functions that already follow block ownership boundaries,
- a documented list of which arrays may be written concurrently and which may not.

This phase is about correctness scaffolding, not speed.

### Phase 4: Add threaded block execution around the shared kernel

Deliverables:

- `parallel/threaded.jl`,
- block chunking and task spawning,
- per-task local accumulators,
- serial reduction of local results,
- one package option to select serial or threaded execution.

This phase should preserve the serial algorithm while changing only execution policy.

### Phase 5: Add regression tests and scaling benchmarks

Deliverables:

- serial versus threaded regression tests,
- test coverage for output structure and restart behavior,
- benchmark scripts for thread scaling,
- baseline comparisons on the example data.

This phase makes later optimization safe.

### Phase 6: Optimize the threaded sampler core

Only after the package structure and tests exist should the main performance work begin.

Targets:

- inner marker-update allocation reduction,
- block scheduling and chunk-size tuning,
- reduction overhead,
- summary accumulation cost,
- restart I/O efficiency,
- GC and memory-pressure tuning when many threads are active.

## Recommendation

The right way to evolve this repository is:

1. keep the current non-MPI package workflow as the single source of truth,
2. refactor that workflow into reusable modules,
3. redesign block data and summaries so they are thread-safe,
4. add Julia multi-threading as a thin execution layer around the shared kernel,
5. validate serial and threaded behavior against each other,
6. then optimize the threaded core.

This path matches the actual shape of the current code, avoids maintaining a second MPI implementation, and gives a practical route from a working research repository to a package-ready scientific software project with a clear shared-memory parallel strategy.