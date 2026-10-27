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
  --gscale_structure correlated \
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
- `--gscale_structure correlated` (`independent` keeps empirical Gscale off-diagonals at zero)
- `--report_pleiotropic_qtl_effect_matrix false`
- `--output_mcmc_delta false`

Argument notes:

- `annot_dict` should match the JLD2 filename stem, for example `anno_matrix_dict` for `anno_matrix_dict.jld2`.
- `starting_value_dir` and `gscale_value_dir` are required by the CLI even for a fresh run; use placeholders such as `XXX` when they are not needed.
- `n_con` is only meaningful for `group_dirichlet`, where the first `n_con` annotation columns after `SNP` are treated as continuous.

## Two-Trait Annotated BayesR

Set

```text
--annotation_prior_model marker_probit_tree_bayesr
```

to use the seven-step annotated two-trait BayesR model. This mode keeps one
latent two-trait marker-effect covariance matrix and assigns each trait a class
in `1:4`, with variance multipliers `(0, 0.01, 0.1, 1.0)`. The first three
probit steps model the `00/10/01/11` activity pattern; four additional steps
model the small/medium/large class within each active trait.

As with `marker_probit_tree`, all annotation columns are prior features,
`estimate_pi` is forced on, and continuation is currently disabled. With
`--output_mcmc_delta true`, `mcmc_Delta1.txt` and `mcmc_Delta2.txt` contain
BayesR class labels `1:4` rather than binary indicators. The run also writes
`annotation_bayesr_joint_prior.txt`, containing posterior summaries of the 16
genome-average annotation-induced joint prior probabilities.

## Optional Gene-Level Genetic Covariance

For a non-MPI `group_dirichlet` run, an exact binary annotation category can
optionally be subdivided using a tab-, comma-, or whitespace-delimited
SNP-to-gene file:

```text
SNP	GENE
rs123	GENE1
```

The `SNP` and `GENE` columns may appear in either order; additional columns are
ignored. Rows with missing or blank SNP or gene identifiers are skipped, and
duplicate SNP–gene pairs are retained once.

Enable the calculation with the required mapping/category pair and the optional
weighting mode:

```bash
--snp-gene-map path/to/snp_gene_map.tsv \
--gene-category exact_annotation_name \
--gene-mapping-weight full
```

`--gene-mapping-weight` accepts `full` (the default) or `fractional`; it can be
omitted when the default is wanted.
Fractional mode multiplies each trait effect by `1/sqrt(k)` when a retained SNP
maps to `k` genes. Gene covariance is calculated from retained `alpha` effects
within each LD block; cross-products are not taken across blocks.

This feature requires an annotation JLD2 created by the updated
`build_annotation_dict`, which stores `my_snp_ids_dict` alongside the unchanged
`my_anno_matrix_dict`. Older annotation JLD2 files continue to work for runs
that do not enable gene covariance.

Enabled runs write
`gene_genetic_covariance_<category>_<weight_mode>.txt` and a matching
`_qc.txt` file in the analysis directory.

The focused implementation checks are in `test/gene_covariance.jl`. A separate
deterministic recovery benchmark in `test/gene_covariance_simulation.jl`
simulates two LD blocks, a gene spanning both blocks, signed effects for two
traits, residual noise, and one SNP mapped to two genes. It runs the actual
non-MPI MCMC workflow in both full and fractional modes and compares posterior
mean gene covariance with the covariance calculated from the simulated effects.
Run only this benchmark with:

```bash
julia --project=. -e 'include("test/gene_covariance_simulation.jl")'
```

The benchmark fixes the marker-effect prior, residual variance, and inclusion
probabilities to isolate gene-covariance recovery. It is an implementation
validation, not a general assessment of model calibration with estimated
hyperparameters.

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
