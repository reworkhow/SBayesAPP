# SBayesAPP

SBayesAPP contains Julia implementations of a bivariate BayesC sampler that uses SNP annotations to model annotation-specific sharing patterns across two traits.

This repository now focuses on a single package-backed execution path that can be optimized further with Julia threads:

- `scripts/run_nonmpi.jl`: primary command-line entrypoint for local and future threaded runs.

The checked-in example under `example/` is set up for this single-process path.

## Repository layout

- `src/`: Julia package source files.
- `scripts/`: Julia workflow entry scripts.
- `script/`: convenience shell launchers.
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
bash script/run.sh
```

If you prefer to call the package entrypoint directly, use:

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

To prepare for thread-based acceleration, run the same entrypoint with Julia threads enabled:

```bash
julia --project=. --threads auto scripts/run_nonmpi.jl \
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

- The workflow creates `analysis_path` automatically if it does not exist.
- `annot_file` is interpreted relative to `data_path`, which is why the example uses `annotation_df.txt` and keeps a copy inside `example/SBayesAPP_input_first10blks/`.
- The optional `--n_con` argument tells SBayesAPP how many annotation columns, starting from the left after `SNP`, should be treated as continuous. Any remaining annotation columns are treated as categorical.
- Non-MPI commands can optionally add `--estimate_vare`, `--estimate_vara`, `--estimate_pi`, `--estimate_gscale`, `--estgscale_iter`, `--report_pleiotropic_qtl_effect_matrix`, and `--output_mcmc_delta`.
- `--gscale_value_dir` is the non-MPI directory used to load fixed `scale_G*.txt` files when continuing with saved Gscale values.
- Set `--report_pleiotropic_qtl_effect_matrix false` to skip the pleiotropic QTL effect-matrix reporting outputs (`MCMC_samples_marker_effects_variance.txt`, `estA*.txt`, `mcmcAtruecor_c.txt`, and related summaries).
- Set `--output_mcmc_delta false` to skip writing `mcmc_Delta*.rank*.txt`; restart output still writes `last_sample_delta/` from the current `deltaArray`.
- If `--estimate_vare false`, `starting_value_dir` should point to a directory containing `estR.txt`.
- `STARTING_VALUE_DIR` and `GSCALE_VALUE_DIR` are placeholders for continuation or fixed-Gscale workflows. They are not used in the provided fresh example run.

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

For a detailed walkthrough of the current non-MPI workflow, see `code_summary.md`.

## Current caveats

- The current non-MPI workflow accepts a seed argument but does not currently call `Random.seed!`.
- Mixed annotation runs are supported when the first `nCon` annotation columns are continuous and the remaining annotation columns are categorical.
- Thread-based speedups still need to be implemented inside the current non-MPI workflow; removing the MPI path is a packaging cleanup rather than a finished threading refactor.

## Next documentation target

The next useful step after this README and the code summary would be to document:

- how the transformed input JLD2 files are generated,
- the exact meaning of each summary output,
- which parts of the sampler are safe to optimize without changing model behavior.