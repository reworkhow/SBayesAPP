# Both Supported Models

This document shows how to run the bundled example data with the two supported annotation prior models:

- `group_dirichlet`
- `marker_probit_tree`

For the `group_dirichlet` quick start, see `readme.md`.

## One-time setup

From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
git lfs pull
```

The example JLD2 files are stored with Git LFS, so `git lfs pull` is required after cloning.

## Shared example paths

These commands use the smaller first-10-block example to keep runtime manageable:

- `DATA_PATH=example/SBayesAPP_input_first10blks/`
- `ST_PATH=example/ST_res/`

Use a separate `ANALYSIS_PATH` for each model so outputs do not overwrite each other.

## Threading recommendation

For multi-threaded runs, it is better to submit `bash_script/run_slurm.sh` than to run `bash_script/run.sh` directly.

`bash_script/run.sh` reads the Julia thread count from `SLURM_CPUS_PER_TASK`, so outside Slurm it usually falls back to one thread.

If you want to run the Julia command directly, use an explicit thread count such as `--threads 3` instead of relying on `--threads auto`.

## Model 1: `group_dirichlet`

You can run `group_dirichlet` either with the bundled shell script or by calling Julia directly.

### Shell script

```bash
bash bash_script/run.sh
```

`bash_script/run.sh` currently follows the `group_dirichlet` path because it does not override the default `--annotation_prior_model` argument.

For a threaded run through the checked-in shell workflow, submit `bash_script/run_slurm.sh` after setting the desired `#SBATCH --cpus-per-task` value.

### Direct Julia command

```bash
julia --project=. --threads 3 scripts/run_nonmpi.jl \
  --data_path example/SBayesAPP_input_first10blks/ \
  --analysis_path example/SBayesAPP_res_first10blks_group_dirichlet/ \
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

`group_dirichlet` uses `--n_con` if you have continuous annotation columns at the left side of the annotation table.

## Model 2: `marker_probit_tree`

Use the same entrypoint, but switch the annotation prior model and write into a different output directory:

```bash
julia --project=. --threads 8 scripts/run_nonmpi.jl \
  --data_path example/SBayesAPP_input_first10blks/ \
  --analysis_path example/SBayesAPP_res_first10blks_marker_probit_tree/ \
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
  --annotation_prior_model marker_probit_tree \
  --is_continue false \
  --estimate_vare true \
  --estimate_vara true \
  --estimate_pi true \
  --estimate_gscale true \
  --estgscale_iter 500 \
  --report_pleiotropic_qtl_effect_matrix false \
  --output_mcmc_delta false
```

Important model-specific behavior:

- `marker_probit_tree` ignores `--n_con`; all annotation columns are treated as prior features.
- `marker_probit_tree` forces `estimate_pi=true`.
- `marker_probit_tree` does not use continuation mode, so keep `--is_continue false`.

For the full argument list and the package CLI defaults, see `readme.md`.

## If you prefer the shell-script workflow for both models

Keep `bash_script/run.sh` for `group_dirichlet`, and run `marker_probit_tree` with the direct Julia command above.

That is the safest option because the checked-in shell script is currently preconfigured for the default `group_dirichlet` path.

## Slurm submission

For either model, you can submit from `bash_script/` once the launcher matches the command you want to run:

```bash
cd bash_script
sbatch run_slurm.sh
```

The current checked-in `run_slurm.sh` simply calls `run.sh`, so if you want a Slurm job for `marker_probit_tree`, create a separate launcher or update the command inside `bash_script/run.sh` first.