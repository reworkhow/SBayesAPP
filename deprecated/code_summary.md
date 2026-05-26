# SBayesAPP Code Summary

## What the code does

SBayesAPP currently exposes a package-backed non-MPI bivariate BayesC-style sampler for transformed summary-statistics inputs. The public entrypoint is `SBayesAPP.run_nonmpi(config)`, and the thin CLI wrapper is `scripts/run_nonmpi.jl`.

The current code supports two annotation prior models:

- `:group_dirichlet`: the original multi-class annotation model. Each annotation category has its own four-state mixture prior and its own marker-effect covariance matrix.
- `:marker_probit_tree`: a marker-specific prior model. The sampler uses a single marker-effect covariance matrix for SNP effects, while annotations feed a three-step probit tree that generates SNP-specific probabilities for the four two-trait inclusion states.

Across both modes, the workflow estimates:

- latent Gaussian marker effects `betaArray`
- inclusion indicators `deltaArray`
- realized marker effects `alphaArray = delta .* beta`
- marker-effect covariance matrices `A_vec`
- genetic covariance summaries `G_vec` and `Gtotal`
- mixture probabilities `Pi`
- block-specific residual covariance matrices `R_blk`, when residual estimation is enabled

## Current module layout

The current non-MPI path is organized around the following files.

- `src/SBayesAPP.jl`: module entrypoint, CLI parsing, command construction, and public `run_nonmpi`.
- `src/config/types.jl`: `NonMPIConfig`, `MarkerProbitTreeState`, and annotation-prior-model validation.
- `src/io/inputs.jl`: annotation metadata loading, block-data loading, and saved effect-state loading.
- `src/io/outputs.jl`: posterior summary writers and restart writers.
- `src/model/priors.jl`: startup `Pi` and `Gprior_vec` construction from single-trait results.
- `src/model/block_setup.jl`: block reindexing, effect-design preparation, and marker-probit-tree design assembly.
- `src/model/initial_state.jl`: parameter initialization and residual-state initialization.
- `src/model/continuation_state.jl`: continuation correction for loaded effects.
- `src/model/mcmc_setup.jl`: posterior-moment allocation and MCMC sample-file setup.
- `src/model/marker_probit_tree.jl`: tree-specific prior state, liability sampling, and prior reconstruction.
- `src/workflows/nonmpi.jl`: the main non-MPI sampler.

## Entry points and configuration

`scripts/run_nonmpi.jl` activates the package, calls `parse_nonmpi_args(ARGS)`, and forwards the resulting `NonMPIConfig` into `run_nonmpi(...)`.

The required named options are:

1. `--data_path`
2. `--analysis_path`
3. `--n_iter`
4. `--annot_file`
5. `--annot_dict`
6. `--out_freq`
7. `--starting_value_dir`
8. `--gscale_value_dir`
9. `--st_path`
10. `--thin`
11. `--n1`
12. `--n2`
13. `--is_continue`

Optional flags include `--seed`, `--n_con`, `--annotation_prior_model`, `--estimate_vare`, `--estimate_vara`, `--estimate_pi`, `--estimate_gscale`, `--estgscale_iter`, `--report_pleiotropic_qtl_effect_matrix`, and `--output_mcmc_delta`.

Important runtime normalization happens before sampling starts.

- `annotation_prior_model` is validated against `(:group_dirichlet, :marker_probit_tree)`.
- In `:marker_probit_tree`, `estimate_pi=false` is ignored and `estimate_pi` is forced on.
- In `:marker_probit_tree`, `is_continue=true` is ignored and the run always starts fresh.
- In `:marker_probit_tree`, `n_con` is ignored and all annotation columns are treated as prior features.
- `estimate_Gscale` is only active when `estimate_vara=true`.

`run_nonmpi(config; dry_run=true)` returns the constructed CLI command instead of launching the workflow.

## Core state representation

For two traits, every SNP belongs to one of four mixture states:

- `[0, 0]`: excluded from both traits
- `[1, 1]`: included in both traits
- `[1, 0]`: included only in trait 1
- `[0, 1]`: included only in trait 2

Those states are stored as dictionaries keyed by length-2 float tuples. The canonical ordering is fixed by `pi_key_order()`:

1. `(0.0, 0.0)`
2. `(1.0, 1.0)`
3. `(1.0, 0.0)`
4. `(0.0, 1.0)`

The main runtime objects are:

- `betaArray`: latent Gaussian marker effects, stored separately for each trait.
- `deltaArray`: binary inclusion indicators with the same layout as `betaArray`.
- `alphaArray`: realized marker effects, computed elementwise as `delta .* beta`.
- `Pi`: a vector of state-probability dictionaries.
- `A_vec`: per-category covariance matrices for latent `beta`.
- `Ainv_vec`: cached inverses of `A_vec`, used in marker updates.
- `R_blk`: per-block residual covariance matrices.
- `G_vec`: per-category genetic covariance summaries from `X * alpha`.

The array layout depends on the annotation prior model.

- In `:group_dirichlet`, `nCategory = nCat + nCon`, and each trait vector in `betaArray`, `deltaArray`, and `alphaArray` has length `my_nsnp * nCategory` in category-major order.
- In `:marker_probit_tree`, `nCategory = 1`, so the effect arrays have length `my_nsnp`. Annotation features no longer create extra effect categories; they only affect prior probabilities.

Tree mode also carries an explicit `MarkerProbitTreeState`, which stores:

- the global annotation design matrix with an intercept
- step-specific regression coefficients and their posterior moments
- step-specific shrinkage variances
- latent liabilities and linear predictors for the three probit steps
- `snp_pi`, the current SNP-specific four-state prior probabilities

In tree mode, `Pi[1]` is only the genome-average summary of `snp_pi` used for output summaries and last-sample files. Marker updates read the marker-specific prior directly from `MarkerProbitTreeState.snp_pi`.

## Annotation prior models

### `:group_dirichlet`

This is the original multi-class annotation formulation.

- Each annotation effect category has its own `Pi[cat]` and `A_vec[cat]`.
- SNPs can participate in multiple effect categories depending on the annotation mask.
- The first `n_con` annotation columns are treated as continuous features. For those columns, the code builds weighted transformed design matrices and forces the sampling mask to `true`.
- Mixture probabilities are updated with a four-state Dirichlet draw from the current per-category state counts.

### `:marker_probit_tree`

This mode keeps the effect model simple and moves annotation dependence into the prior.

- The sampler uses a single effect category for all SNPs.
- The annotation matrix is turned into one global design matrix with an intercept.
- A three-step probit tree defines marker-specific state probabilities.

The three steps are:

1. `step1_zero_vs_active`: zero state versus any active state
2. `step2_11_vs_singleton`: shared effect `(1, 1)` versus singleton states among active SNPs
3. `step3_10_vs_01`: trait-1-only versus trait-2-only among singleton states

Let `Phi` denote the standard normal CDF. For each SNP, the three step-specific linear predictors `mu1`, `mu2`, and `mu3` are converted into branch probabilities `p1 = Phi(mu1)`, `p2 = Phi(mu2)`, and `p3 = Phi(mu3)`. The four-state prior is then obtained by following the tree:

- `P(0,0) = 1 - p1`: stop at step 1 and stay inactive
- `P(1,1) = p1 * p2`: pass step 1, then choose the shared-effect branch at step 2
- `P(1,0) = p1 * (1 - p2) * p3`: pass step 1, go to the singleton branch at step 2, then choose trait 1 at step 3
- `P(0,1) = p1 * (1 - p2) * (1 - p3)`: pass step 1, go to the singleton branch at step 2, then choose trait 2 at step 3

Step coefficients are sampled from Gaussian conditionals after drawing truncated-normal liabilities. The intercept is unpenalized. If the tree design matrix includes annotation columns beyond the intercept, those annotation coefficients are shrunk with step-specific variances, and the variances are resampled each iteration.

## Input contract

### Shared block-wise inputs

The non-MPI workflow expects transformed block-wise inputs prepared upstream.

- `TransformedX_dict.jld2`: block id to transformed genotype basis `X`
- `TransformedY_dict.jld2`: block id to transformed response vectors for the two traits
- `blkSNPsIndex_dict.jld2`: per-block SNP indices
- `blkIDs.txt`: block ids to analyze
- `nGWAS_dict.jld2`: per-block sample-size information
- `annot_dict.jld2`: per-block annotation matrix dictionary
- `annot_file`: annotation table used to derive annotation names, sizes, and types

The example directory `example/SBayesAPP_input_first10blks/` matches the current loader contract.

### Single-trait startup inputs

The workflow seeds several priors from a single-trait result directory `st_path`.

- `Trait1/mean_pi.txt`
- `Trait2/mean_pi.txt`
- `Trait1/mean_varg_total.txt`
- `Trait2/mean_varg_total.txt`

These files are used to build the startup four-state prior `startPi` and the initial marker-effect covariance prior scale `Gprior_vec`.

### Annotation handling by model

The same raw annotation files are interpreted differently by the two models.

- `:group_dirichlet`: `load_annotation_metadata(...; nCon=n_con)` treats the first `n_con` columns as continuous and the rest as categorical. `prepare_block_state!` builds weighted transformed designs and a Boolean sampling mask derived from the annotation matrix.
- `:marker_probit_tree`: `n_con` is ignored. `prepare_marker_probit_tree_block_state!` builds a single all-marker sampling mask, uses the unweighted transformed design, and assembles a global annotation design matrix with an intercept. Annotation values affect prior probabilities only.

Continuation currently belongs only to the group-dirichlet path. When continuation is active, the code reloads `betaArray`, `deltaArray`, `Pi`, `A_vec`, and `R_blk`, reconstructs `alphaArray`, and subtracts those effects back out of the transformed responses before resuming.

## Setup phase

Before the MCMC loop starts, the current workflow does the following.

1. `build_nonmpi_run_context` loads annotation metadata and normalizes model-specific settings.
2. `build_start_pi` constructs the startup four-state prior from the single-trait `mean_pi` files.
3. In `:group_dirichlet`, `build_gprior_vec` is computed immediately from annotation sizes. In `:marker_probit_tree`, the `Gprior_vec` construction is deferred until `my_nsnp` is known, then built for the single effect category.
4. `load_nonmpi_block_data` loads transformed inputs, block ids, block SNP indices, annotation matrices, and sample-size dictionaries.
5. The workflow chooses the sampler dimensions:
   - `nCategory = annotation_metadata.nCat + nCon` for `:group_dirichlet`
   - `nCategory = 1` and `nLoci_annot = [my_nsnp]` for `:marker_probit_tree`
6. `initialize_nonmpi_parameter_state` initializes `A_vec`, `Ainv_vec`, `Pi`, `Rprior`, `scale_G_vec`, and related hyperparameters.
7. `prepare_block_state!` or `prepare_marker_probit_tree_block_state!` reorders SNP indices across blocks and builds the design objects needed inside the hot loop.
8. If tree mode is active, `initialize_marker_probit_tree_state` creates the intercept-plus-annotation tree state and initializes `snp_pi` from `startPi`.
9. The workflow initializes `betaArray`, `alphaArray`, `deltaArray`, posterior mean buffers, optional delta storage, and the MCMC sample files.

## Main MCMC loop

The main iteration loop is `for iter = 1:nIter` inside `run_nonmpi_sampler!`.

### Shared marker update

For each iteration, block, marker, and active category, the sampler:

1. Builds the working statistics `w = x' * (y_corr + x * alpha_old)` for both traits.
2. Samples the inclusion indicator and latent effect for each trait conditional on the current other trait state.
3. Updates the working residual vector in place.
4. Writes the new `beta`, `delta`, and `alpha` values back into the global arrays.
5. Increments the four-state counts used for the later `Pi` update.

The prior term inside the delta update comes through `log_marker_state_prior(...)`.

- In `:group_dirichlet`, it reads `Pi[cat]`.
- In `:marker_probit_tree`, it reads the marker-specific row from `MarkerProbitTreeState.snp_pi`.

This nested block-marker-category-trait loop remains the dominant compute cost in the current implementation.

### Shared variance and covariance updates

After each block update, the workflow:

- reconstructs block-level total genetic variance from `X * alpha` whenever either `do_thin` is true for that iteration or `estimate_vare=true`
- records category-level genetic variance and covariance summaries only on thinning iterations
- records per-category `alpha'alpha` sums only when residual variance estimation is enabled
- optionally samples a block-specific residual covariance `R_blk[b]`

If residual estimation is enabled, each block first gets an inverse-Wishart residual draw. The code then applies a diagonal-threshold heuristic trait by trait: it computes `sum(ssq_blk_cat[b][traiti, :]) / totalvarg_blk[b][traiti, traiti]`, keeps the sampled diagonal only when that ratio exceeds `1.1`, and otherwise resets the diagonal to `1 / nInd[traiti]`. After that, it rebuilds the off-diagonal entry from the sampled residual correlation and the retained diagonals. The current code averages those block-specific matrices only when recording posterior summaries; it does not replace all `R_blk[b]` with a common matrix for the next iteration.

When thinning is active, the workflow updates running means for:

- `meanAlpha`
- `meanA`, when pleiotropic marker-effect covariance reporting is enabled
- `meanB`
- `meanG`
- `meanGtotal`
- `meanR`, when residual estimation is enabled

### Model-specific `Pi` update

The two prior models diverge at the mixture-probability update step.

- `:group_dirichlet`: each annotation category samples a new four-state Dirichlet draw from `nLoci_array_vec[cat] .+ 1`, and those values overwrite `Pi[cat]`.
- `:marker_probit_tree`: `update_marker_probit_tree_priors!` converts the current `deltaArray` into three binary response vectors, samples liabilities for each step, updates the step coefficients and step variances, rebuilds SNP-specific mixture probabilities, and stores the genome-average mixture probabilities back into `Pi[1]` for summaries.

### Shared covariance-hyperparameter update

The code recomputes `SSE_vec = beta'beta` per category on every iteration.

- During the first `estGscale_iter` iterations, `scale_G_vec` can be estimated from running `SSE` averages if `estimate_Gscale=true`.
- If `estimate_vara=true`, each `A_vec[cat]` is sampled from its inverse-Wishart full conditional posterior, using `scale_G_vec[cat] + SSE_vec[cat]` as the scale matrix for the full conditional posterior distribution, and `Ainv_vec` is refreshed for the next MCMC iteration's marker updates.

### Checkpointing and sample writes

After burn-in, the workflow writes different outputs on two schedules.

- At every `out_freq` checkpoint, it appends the current `Pi` samples and `A_vec` samples, and optionally saves `delta` trajectories.
- At every thinning point, it appends total and category-level genetic covariance samples, and optionally the pleiotropic realized-effect covariance samples.
- At the last scheduled save iteration, it writes the current samples to the last-sample files for `betaArray`, `R_blk`, `A_vec`, `Pi`, and `deltaArray`.

## Final outputs

### Common outputs

The current workflow writes the following summary files when the relevant options are enabled.

- `annotationName.txt`
- `MCMC_samples_pi.txt`
- `MCMC_samples_beta_effects_variance.txt`
- `MCMC_samples_genetic_effects_variance.txt`
- `MCMC_samples_total_genetic_effects_variance.txt`
- `estB*.txt` and `estB_std*.txt`
- `estG*.txt` and `estG_std*.txt`
- `estGtotal.txt` and `estGtotal_std.txt`
- `estGcor_total.txt` and `estGcor_total_std.txt`
- `estR.txt` and `estR_std.txt`, when `estimate_vare=true`
- `estPi*.txt` and `estPi_std*.txt`, when `estimate_pi=true`
- `meanAlpha*.rank*.txt`
- `mcmc_Delta*.rank*.txt`, when `output_mcmc_delta=true`
- last-sample files under `last_sample_R_blk/`, `beta_effect_var_matrices_last_sample/`, `pi_last_sample/`, and `last_sample_delta/`

If `report_pleiotropic_qtl_effect_matrix=true`, the workflow also writes:

- `MCMC_samples_marker_effects_variance.txt`
- `estA*.txt` and `estA_std*.txt`

### `:group_dirichlet`-specific outputs

Additional covariance and correlation trace files are only written in group-dirichlet mode. This includes the total MCMC trace files `mcmcGcov_total.txt` and `mcmcGcor_total.txt`, as well as the category-level correlation files below. The total posterior summary files listed above, such as `estGcor_total.txt`, are still written in both models.

- `mcmcGcov_total.txt`
- `mcmcGcor_total.txt`
- `mcmcGcov_c.txt`
- `mcmcGcor_c.txt`
- `mcmcAtruecor_c.txt`, when pleiotropic reporting is enabled
- `estBcor.txt` and `estBcor_std.txt`
- `estAcor.txt` and `estAcor_std.txt`, when pleiotropic reporting is enabled
- `estGcor.txt` and `estGcor_std.txt`

### `:marker_probit_tree`-specific outputs

Tree mode suppresses the category-level correlation summary files as above, and instead writes the probit regression summaries:

- `annotation_probit_coefficients.txt`
- `annotation_probit_coefficients_std.txt`

In this mode, `estPi1.txt` is the posterior mean of the genome-average of marker-specific mixture probabilities, not a true category-level mixture prior shared identically by all markers.

## Important implementation details

- `run_nonmpi_workflow` calls `Random.seed!(seed)` before building the sampler context. If the `seed` argument is omitted, the code uses the default seed value `123`.
- `analysis_path` is created before outputs are written.
- `Pi` output ordering is stabilized by `pi_key`, `pi_key_order`, and `write_pi_dict`, not by raw dictionary iteration order.
- Tree-mode startup requires a valid `startPi` with positive shared-state mass and positive active mass for both traits.
- If `gscale_value_dir` is supplied on continuation, `scale_G_vec` is loaded from disk and no longer re-estimated.

## Optimization-relevant observations

These points matter if the next step is code or repository optimization.

- The dominant cost is still the nested loop over blocks, markers, categories, and traits.
- `:marker_probit_tree` reduces the category dimension to one, but it adds an extra dense annotation-regression update over all markers each iteration.
- The code still allocates temporary vectors and slices inside hot loops, especially around residual updates, per-marker state reconstruction, and category slices.
- Restart and summary output is still text-heavy, which is convenient for inspection but expensive for large analyses.
- Memory usage still scales with `my_nsnp * nCategory` per trait for `betaArray`, `deltaArray`, and `alphaArray`, so the group-dirichlet mode remains sensitive to the number of effect categories.

That makes this documentation useful groundwork for later refactoring: most performance work should preserve the current update order while reducing allocation, duplication, and unnecessary per-iteration IO.