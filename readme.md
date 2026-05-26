# SBayesAPP Example:

This document shows the shortest path to run SBayesAPP with the bundled example data using the `group_dirichlet` annotation prior model.

If you also want to run `marker_probit_tree`, see `deprecated/readme_models.md`.

## Before you run

From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
git lfs pull
```

`git lfs pull` is required because the example JLD2 inputs in `example/` are stored with Git LFS.

## Example data used by the default launcher

The default shell launcher in `bash_script/run.sh` uses the smaller first-10-block example so the run finishes faster:

- `example/SBayesAPP_input_first10blks/`: transformed input dictionaries.
- `example/ST_res/`: single-trait summaries used to initialize the bivariate run.
- `example/SBayesAPP_res_first10blks/`: output directory created by the launcher.

The shell script currently uses the default CLI model, which is `group_dirichlet`.

## Quick start with the bundled shell script

From the repository root:

```bash
bash bash_script/run.sh
```

This script currently runs with:

- `DATA_PATH=example/SBayesAPP_input_first10blks/`
- `ANALYSIS_PATH=example/SBayesAPP_res_first10blks/`
- `NITER=500`
- `ST_PATH=example/ST_res/`
- Julia threads set from `SLURM_CPUS_PER_TASK` when available, otherwise `1`

If you want to change the input path, output path, iteration count, or sample sizes, edit the variables at the top of `bash_script/run.sh`.

If you run `bash bash_script/run.sh` directly in a normal shell, it will usually run with a single thread because `SLURM_CPUS_PER_TASK` is not set.

## Recommended multi-threaded run

If you want a multi-threaded example run, the recommended path is:

```bash
cd bash_script
sbatch run_slurm.sh
```

`bash_script/run_slurm.sh` is the better entrypoint for multi-threading because `bash_script/run.sh` reads the Julia thread count from `SLURM_CPUS_PER_TASK`.

To change the thread count, edit `#SBATCH --cpus-per-task` in `bash_script/run_slurm.sh` before submitting.

## Equivalent Julia command

If you want to run the Julia entrypoint directly, set an explicit thread count instead of relying on `auto`. For example:

```bash
julia --project=. --threads 3 scripts/run_nonmpi.jl \
  --data_path example/SBayesAPP_input_first10blks/ \
  --analysis_path example/SBayesAPP_res_first10blks/ \
  --n_iter 500 \
  --seed 42 \
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
  --annotation_prior_model group_dirichlet \
  --is_continue false \
  --estimate_vare true \
  --estimate_vara true \
  --estimate_pi true \
  --estimate_gscale true \
  --estgscale_iter 500 \
  --report_pleiotropic_qtl_effect_matrix false \
  --output_mcmc_delta false
```

Replace `3` with the number of threads you want to use on your machine.

`annot_file` is resolved relative to `data_path`, so `annotation_df.txt` must exist inside `example/SBayesAPP_input_first10blks/` for this command.

## Run through Slurm

The Slurm wrapper expects to be submitted from the `bash_script/` directory so it can find `run.sh`.

```bash
cd bash_script
sbatch run_slurm.sh
```

Before submitting, adjust resource settings in `bash_script/run_slurm.sh` if needed.

## Required arguments and CLI defaults

The SBayesAPP command requires these arguments every time:

- `--data_path`
- `--analysis_path`
- `--n_iter`
- `--annot_file`
- `--annot_dict`
- `--out_freq`
- `--starting_value_dir`
- `--gscale_value_dir`
- `--st_path`
- `--thin`
- `--n1`
- `--n2`
- `--is_continue`

The package-level CLI defaults are:

- `--seed 123`
- `--n_con 0`
- `--annotation_prior_model group_dirichlet`
- `--estimate_vare true`
- `--estimate_vara true`
- `--estimate_pi true`
- `--estimate_gscale true`
- `--estgscale_iter 500`
- `--report_pleiotropic_qtl_effect_matrix false`
- `--output_mcmc_delta false`

The bundled `bash_script/run.sh` overrides some of those defaults for the example run, including `seed=42` and `n_iter=500`.

Argument notes:

- `annot_file` is interpreted relative to `data_path`.
- `annot_dict` should match the JLD2 filename stem, for example `anno_matrix_dict` for `anno_matrix_dict.jld2`.
- `starting_value_dir` and `gscale_value_dir` are still required by the CLI even for a fresh run; use placeholders such as `XXX` when they are not needed.
- `n_con` is only meaningful for `group_dirichlet`, where the first `n_con` annotation columns after `SNP` are treated as continuous.

## Required files for this example

The example data directory must contain:

- `TransformedX_dict.jld2`
- `TransformedY_dict.jld2`
- `blkSNPsIndex_dict.jld2`
- `blkIDs.txt`
- `nGWAS_dict.jld2`
- `anno_matrix_dict.jld2`
- `annotation_df.txt`

The single-trait initialization directory must contain:

- `Trait1/mean_pi.txt`
- `Trait2/mean_pi.txt`
- `Trait1/mean_varg_total.txt`
- `Trait2/mean_varg_total.txt`

## Main outputs

The example run writes results under `example/SBayesAPP_res_first10blks/`, including:

- MCMC traces such as `MCMC_samples_pi.txt`
- posterior summaries such as `estPi*.txt`, `estG*.txt`, `estB*.txt`, and `estGtotal.txt`
- last mcmc sample files such as `last_mcmc_betaArray*.rank*.txt`, `pi_last_sample/`, `beta_effect_var_matrices_last_sample/`, and `last_sample_delta/`

For both supported annotation-prior models, see `deprecated/readme_models.md`.