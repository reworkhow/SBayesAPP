# SBayesAPP

SBayesAPP contains Julia implementations of a bivariate BayesC sampler that uses SNP annotations to model annotation-specific sharing patterns across two traits.

This repository currently has two main execution paths:

- `src/app_nonMPI.jl`: single-process sampler for a local or small example run.
- `src/app_MPI.jl`: MPI-enabled sampler for distributed runs where LD blocks are split across ranks.

The checked-in example under `example/` is set up for the non-MPI path.

## Repository layout

- `src/`: Julia source files and package entrypoints.
- `scripts/`: Julia workflow entry scripts.
- `script/`: compatibility shell launcher.
- `example/`: demonstration inputs and outputs.
- `context.md`: short repository note describing the intended example workflow.

## Environment setup

From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This repository now includes a `Project.toml`, so dependency resolution can happen through Julia directly.

## Included example

The example directory includes:

- `example/SBayesAPP_input_first10blks/`: transformed input for the first 10 LD blocks.
- `example/ST_res/`: single-trait outputs used to initialize priors.
- `example/annotation_df.txt`: annotation table.
- `example/SBayesAPP_res_first10blks/`: example output directory.

The current shell launcher is intended to run the single-process example and write results into the example output directory.

## Quick start

From the repository root:

```bash
julia --project=. scripts/run_example.jl
```

If you prefer the compatibility shell wrapper, this still works:

```bash
bash script/run.sh
```

## Direct non-MPI run

If you want to call the current entrypoint directly, the equivalent command is:

```bash
julia --project=. src/app_nonMPI.jl \
  --data_path example/SBayesAPP_input_first10blks/ \
  --analysis_path example/SBayesAPP_res_first10blks/ \
  --n_iter 1000 \
  --seed 42 \
  --nrank 1 \
  --annot_file annotation_df.txt \
  --annot_dict anno_matrix_dict \
  --out_freq 100 \
  --starting_value_dir XXX \
  --gscale_value_dir XXX \
  --st_path example/ST_res/ \
  --thin 50 \
  --n1 300000 \
  --n2 300000 \
  --n_con 0 \
  --is_continue false
```

There is also a package-backed wrapper script that forwards those same named arguments:

```bash
julia --project=. scripts/run_nonmpi.jl \
  --data_path example/SBayesAPP_input_first10blks/ \
  --analysis_path example/SBayesAPP_res_first10blks/ \
  --n_iter 1000 \
  --seed 42 \
  --nrank 1 \
  --annot_file annotation_df.txt \
  --annot_dict anno_matrix_dict \
  --out_freq 100 \
  --starting_value_dir XXX \
  --gscale_value_dir XXX \
  --st_path example/ST_res/ \
  --thin 50 \
  --n1 300000 \
  --n2 300000 \
  --n_con 0 \
  --is_continue false
```

Notes:

- The output directory should exist before launching `src/app_nonMPI.jl` directly.
- `annot_file` is interpreted relative to `data_path`, which is why the example uses `annotation_df.txt` and keeps a copy inside `example/SBayesAPP_input_first10blks/`.
- The optional `--n_con` argument tells SBayesAPP how many annotation columns, starting from the left after `SNP`, should be treated as continuous. Any remaining annotation columns are treated as categorical.
- Non-MPI commands can optionally add `--estimate_vare`, `--estimate_vara`, `--estimate_pi`, `--estimate_gscale`, `--estgscale_iter`, `--report_pleiotropic_qtl_effect_matrix`, and `--output_mcmc_delta`.
- `--gscale_value_dir` is the non-MPI directory used to load fixed `scale_G*.txt` files when continuing with saved Gscale values.
- Set `--report_pleiotropic_qtl_effect_matrix false` to skip the pleiotropic QTL effect-matrix reporting outputs (`MCMC_samples_marker_effects_variance.txt`, `estA*.txt`, `mcmcAtruecor_c.txt`, and related summaries).
- Set `--output_mcmc_delta false` to skip writing `mcmc_Delta*.rank*.txt`; restart output still writes `last_sample_delta/` from the current `deltaArray`.
- If `--estimate_vare false`, `starting_value_dir` should point to a directory containing `estR.txt`.
- `STARTING_VALUE_DIR` and `GSCALE_VALUE_DIR` are placeholders for continuation or fixed-Gscale workflows. They are not used in the provided fresh example run.

## MPI run

The MPI entrypoint is `src/app_MPI.jl`. It expects rank-partitioned input files and a few extra arguments beyond the non-MPI version, including average sample sizes and flags controlling whether hyperparameters are estimated or fixed.

A typical invocation pattern is:

```bash
mpiexec -n <ranks> julia --project=. scripts/run_mpi.jl \
  --data_path <data_path>/ \
  --analysis_path <analysis_path>/ \
  --n_iter <nIter> \
  --seed <seed> \
  --nrank <ranks> \
  --annot_file <annot_file> \
  --annot_dict <annot_dict> \
  --out_freq <outFreq> \
  --starting_value_dir <starting_value_dir> \
  --secondary_starting_value_dir <secondary_starting_value_dir> \
  --st_path <ST_path>/ \
  --thin <thin> \
  --n1 <N1> \
  --n2 <N2> \
  --n_con [nCon] \
  --estimate_pi <estimate_pi> \
  --fixed_hyperparameters [fixed_hyperparameters] \
  --is_continue [is_continue] \
  --chr [chr]
```

This path is intended for pre-sharded data, not for the small example in `example/`.

## Required inputs

For a non-MPI run, the data directory is expected to contain:

- `TransformedX_dict.jld2`
- `TransformedY_dict.jld2`
- `blkSNPsIndex_dict.jld2`
- `blkIDs.txt`
- `nGWAS_dict.jld2`
- `<annot_dict>.jld2`
- `annotation_df.txt` or another annotation file named by the CLI argument

The single-trait initialization directory is expected to contain:

- `Trait1/mean_pi.txt`
- `Trait2/mean_pi.txt`
- `Trait1/mean_varg_total.txt`
- `Trait2/mean_varg_total.txt`

## Output overview

The samplers write three kinds of outputs.

- Running MCMC samples such as `MCMC_samples_pi.txt` and `MCMC_samples_genetic_effects_variance.txt`.
- Posterior summaries such as `estPi*.txt`, `estG*.txt`, `estGtotal.txt`, and `estR.txt`.
- Restart files such as `last_mcmc_betaArray*.rank*.txt`, `pi_last_sample/`, `beta_effect_var_matrices_last_sample/`, and `last_sample_delta/`.

For a detailed walkthrough of how those files are produced and how the two Julia entrypoints differ, see `code_summary.md`.

## Current caveats

- The repository does not yet include a Julia project file for reproducible dependency installation.
- `src/app_nonMPI.jl` parses a seed argument but does not currently call `Random.seed!`.
- Mixed annotation runs are supported when the first `nCon` annotation columns are continuous and the remaining annotation columns are categorical.
- The MPI code contains additional workflow modes, including split-chromosome and fixed-hyperparameter runs, but the checked-in example covers only the non-MPI path.

## Next documentation target

The next useful step after this README and the code summary would be to document:

- how the transformed input JLD2 files are generated,
- the exact meaning of each summary output,
- which parts of the sampler are safe to optimize without changing model behavior.