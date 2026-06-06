# SBayesAPP

SBayesAPP (Summary-data-based Bayesian method leveraging biological Annotations to quantify Pleiotropy and Polygenicity) is a Bayesian framework for dissecting shared genetic architecture between complex traits using GWAS summary statistics and functional annotations.

The method estimates:

- annotation-stratified genetic covariance matrices
- annotation-stratified SNP effect covariance matrices
- annotation-stratified polygenic proportions

SBayesAPP helps distinguish whether coheritability enrichment is driven by many shared variants with weak effects or by fewer variants with stronger pleiotropic effects.

## Requirements

SBayesAPP is written in Julia and uses R scripts for some preprocessing, downstream analysis, and visualization.

Typical requirements:

- Julia 1.11.x
- R 4.1 or newer for analysis/plotting scripts
- Git LFS for bundled example data

From the repository root, instantiate the Julia environment before running examples:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
git lfs pull
```

`git lfs pull` is required because the example JLD2 input files under `example/` are stored with Git LFS.

## Data

The T2D and fasting glucose cell-type analysis data are available here:

https://drive.google.com/drive/folders/1nR_wAj9Hwk1LCrRcFVEr9J-qXiNeaaRp

## Repository Structure

```text
SBayesAPP/
├── Analysis_codes/             # real-data and simulation analysis code used in downstream analyses
├── bash_script/                # shell launchers for local and Slurm example runs
├── deprecated/                 # older notes, model docs, and historical scripts
├── example/                    # bundled example input and output files
├── scripts/                    # Julia command-line entrypoints
├── src/                        # SBayesAPP Julia package source
├── test/                       # package tests
├── Project.toml
├── Manifest.toml
└── README.md
```

## Quick Start

The default shell launcher runs the smaller first-10-block example so it finishes faster.

From the repository root:

```bash
bash bash_script/run.sh
```

The default launcher uses:

- `DATA_PATH=example/SBayesAPP_input_first10blks/`
- `ANALYSIS_PATH=example/SBayesAPP_res_first10blks/`
- `NITER=500`
- `ST_PATH=example/ST_res/`
- Julia threads from `SLURM_CPUS_PER_TASK` when available, otherwise `1`

To change the input path, output path, iteration count, or sample sizes, edit the variables near the top of `bash_script/run.sh`.

## Multi-threaded Slurm Example

For a multi-threaded example run, submit the Slurm wrapper from the `bash_script/` directory:

```bash
cd bash_script
sbatch run_slurm.sh
```

`bash_script/run_slurm.sh` is the preferred entrypoint for multi-threaded example runs because `bash_script/run.sh` reads the Julia thread count from `SLURM_CPUS_PER_TASK`.

To change the thread count, edit `#SBATCH --cpus-per-task` in `bash_script/run_slurm.sh` before submitting.

## Equivalent Julia Command

You can also run the Julia entrypoint directly. Set an explicit thread count instead of relying on `auto`:

```bash
julia --project=. --threads 3 scripts/run_nonmpi.jl \
  --data_path example/SBayesAPP_input_first10blks/ \
  --analysis_path example/SBayesAPP_res_first10blks/ \
  --n_iter 500 \
  --burnin 200 \
  --seed 42 \
  --annot_file annotation_df.txt \
  --annot_dict anno_matrix_dict \
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

Replace `3` with the number of threads you want to use.

`annot_file` is resolved relative to `data_path`, so `annotation_df.txt` must exist inside `example/SBayesAPP_input_first10blks/` for this command.

`--burnin` is optional. If omitted, fresh runs default to trimming 40% of the chain, while continuation runs default to `0`.

## Required Arguments

The SBayesAPP command requires these arguments:

- `--data_path`
- `--analysis_path`
- `--n_iter`
- `--annot_file`
- `--annot_dict`
- `--starting_value_dir`
- `--gscale_value_dir`
- `--st_path`
- `--thin`
- `--n1`
- `--n2`
- `--is_continue`

`--thin` controls both posterior-mean thinning and periodic MCMC sample/checkpoint writes.

Package-level CLI defaults include:

- `--burnin floor(0.4 * n_iter)` for fresh runs, `0` for continuation runs
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

Argument notes:

- `annot_dict` should match the JLD2 filename stem, for example `anno_matrix_dict` for `anno_matrix_dict.jld2`.
- `starting_value_dir` and `gscale_value_dir` are required by the CLI even for a fresh run; use placeholders such as `XXX` when they are not needed.
- `n_con` is only meaningful for `group_dirichlet`, where the first `n_con` annotation columns after `SNP` are treated as continuous.

## Required Example Files

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

## Main Outputs

The example run writes results under `example/SBayesAPP_res_first10blks/`, including:

- MCMC traces such as `MCMC_samples_pi.txt`
- posterior summaries such as `estPi*.txt`, `estG*.txt`, `estB*.txt`, and `estGtotal.txt`
- last-sample files such as `last_mcmc_betaArray*.rank*.txt`, `pi_last_sample/`, `beta_effect_var_matrices_last_sample/`, and `last_sample_delta/`

## Additional Documentation

For notes about supported annotation-prior models and older workflows, see `deprecated/readme_models.md`.
