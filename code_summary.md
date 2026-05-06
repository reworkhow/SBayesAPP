# SBayesAPP Code Summary

## What the code does

SBayesAPP implements a bivariate BayesC-style sampler that uses annotation-defined SNP groups to estimate:

- per-annotation mixture proportions for the four inclusion states across two traits,
- per-annotation SNP effect covariance matrices,
- per-annotation genetic covariance matrices,
- total genetic covariance matrices across the genome,
- residual covariance matrices, when that component is enabled.

The repository now focuses on the package-backed non-MPI workflow, exposed through `scripts/run_nonmpi.jl`.

## Shared model objects

The current non-MPI workflow organizes the sampler around the same parameterization used throughout the package.

- `betaArray`: latent Gaussian marker effects, stored separately for each trait, and within each trait vector the effects are ordered by annotation category, for example `[SNP1_c1 SNP2_c1 SNP1_c2 SNP2_c2]`, where `SNPi_cj` is the effect of marker `i` in category `j`.
- `deltaArray`: binary inclusion indicators for each trait at each SNP ordered by category. The ordering matches `betaArray`.
- `alphaArray`: realized marker effects used in the likelihood, defined as elementwise `delta .* beta`.
- `Pi`: annotation-specific probabilities over the four bivariate inclusion states.
- `A_vec`: annotation-specific covariance matrices for `beta`.
- `R_blk`: block-level residual covariance matrices.
- `G_vec`: annotation-specific genetic covariance estimated from `X * alpha`.

For two traits, each SNP in each category belongs to one of four mixture states:

- `[0, 0]`: excluded from both traits
- `[1, 1]`: included in both traits
- `[1, 0]`: included only in trait 1
- `[0, 1]`: included only in trait 2

The code stores those probabilities in dictionaries keyed by the corresponding length-2 tuples.

## Input contract

### Shared data inputs

The non-MPI workflow expects transformed block-wise inputs that have already been prepared upstream.

- `TransformedX_dict.jld2`: dictionary from block id to transformed genotype basis `X`.
- `TransformedY_dict.jld2`: dictionary from block id to transformed response vectors for trait 1 and trait 2.
- `blkSNPsIndex_dict.jld2`: per-block SNP indices.
- `blkIDs.txt`: block ids to analyze.
- `nGWAS_dict.jld2`: per-block sample size information.
- annotation matrix dictionary: a per-block annotation membership matrix saved as JLD2.
- annotation table: used to determine the annotation order and the number of SNPs in each annotation.

The current example input directory under `example/SBayesAPP_input_first10blks/` matches the current non-MPI loader.

### Single-trait initialization inputs

The workflow sets up several priors from a single-trait result directory `ST_path`, specifically:

- `Trait1/mean_pi.txt`
- `Trait2/mean_pi.txt`
- `Trait1/mean_varg_total.txt`
- `Trait2/mean_varg_total.txt`

These files are used to build:

- the initial mixture prior `startPi`,
- the initial scale of annotation-specific marker effect covariance matrix priors `Gprior_vec`.

### Annotation handling

The current package path supports both categorical-only annotations and mixed continuous-plus-categorical annotations.

- `--n_con` controls how many annotation columns, starting from the left after `SNP`, are treated as continuous.
- Any remaining annotation columns are treated as categorical.
- Annotation metadata is loaded once and then used to build the per-block masks and transformed design views consumed by the sampler.

## `scripts/run_nonmpi.jl`

### Command-line interface

`scripts/run_nonmpi.jl` forwards named CLI arguments into `parse_nonmpi_args(ARGS)` and then calls `run_nonmpi(...)`.

The required options are:

1. `--data_path`
2. `--analysis_path`
3. `--n_iter`
4. `--seed`
5. `--nrank`
6. `--annot_file`
7. `--annot_dict`
8. `--out_freq`
9. `--starting_value_dir`
10. `--gscale_value_dir`
11. `--st_path`
12. `--thin`
13. `--n1`
14. `--n2`
15. `--is_continue`

Optional flags extend that contract with `--n_con`, `--estimate_vare`, `--estimate_vara`, `--estimate_pi`, `--estimate_gscale`, `--estgscale_iter`, `--report_pleiotropic_qtl_effect_matrix`, and `--output_mcmc_delta`.

In practice, `nrank` is expected to be `1` in this file.

### Setup phase

The non-MPI workflow does the following before entering MCMC:

1. Reads the annotation table and derives the ordered annotation names and annotation sizes.
2. Reads single-trait outputs to initialize `startPi` and `Gprior_vec`.
3. Loads per-block transformed data dictionaries from `data_path`.
4. Reindexes SNP positions so that SNPs are sequential across blocks rather than local within each block.
5. Precomputes:
   - `xpx_dict`: per-marker quadratic terms `x'x` for each block.
   - `xArray_dict`: per-block transformed design matrices.
6. Converts the annotation matrix into a Boolean sampling mask.
7. Initializes `betaArray`, `alphaArray`, `deltaArray`, and summary accumulators.

If `is_continue == true`, the code also restores previous samples for:

- `A_vec`
- `Pi`
- `betaArray`
- `deltaArray`
- block-specific `R_blk`

It then subtracts the previously loaded effects from `my_TransformedY_dict` so the residual state is consistent with continuing the chain.

### Main MCMC loop

The heart of the program is `run_nonmpi_sampler!`, whose iteration loop starts at `for iter = 1:nIter`.

Each iteration does the following.

#### 1. Reset per-iteration accumulators

The code zeroes block-level and category-level accumulators for:

- per-block per-category genetic variance and covariance (`varg_blk_cat` and `varg_cov_blk_cat`),
- per-category genetic covariance matrices (`G_vec`),
- per-block total genetic covariance matrices,
- per-category distributions of loci in the four mixture states,
- per-category sum of squares of latent Gaussian marker effects (`SSE_vec`),
- per-block per-category sums of squares of realized marker effects, which are only used to determine the strategy of residual covariance matrix sampling.

#### 2. Loop over LD blocks

For each block:

1. Read the transformed design matrix, working response vectors, annotation mask, and sample sizes.
2. Build the block residual precision `Rinv` from the current block-specific residual covariance.
3. Randomize marker order.

#### 3. Loop over markers and annotation categories

For each marker and each annotation category that includes that marker:

1. Recover the current two-trait state from `betaArray`, `alphaArray`, and `deltaArray`.
2. Compute the scalar working statistic `w = x' (y_corr + x * alpha_old)` for each trait.
3. Randomize trait order.
4. For each trait, compare the posterior support for `delta = 0` versus `delta = 1` using the current sampled values of `Pi[cat]` and covariance precision `Ainv_vec[cat]`.
5. Sample:
   - the inclusion indicator `delta[k]`,
   - the latent Gaussian effect `beta[k]`,
   - the realized effect `alpha[k]`.
6. Update the working residual vector `wArray` in place.
7. Increment the annotation-specific mixture-state counts used later for updating `Pi`.
8. Write the new values back to the global effect arrays.

This is the expensive part of the algorithm. Its computational cost scales with the product of:

- number of iterations,
- number of blocks,
- number of markers per block,
- number of annotation groups,
- number of traits.

#### 4. Compute per-block genetic variance summaries

After updating all markers in a block, the code reconstructs annotation-specific genetic contributions using `X * alpha`.

For each category it first computes block-level:

- per-trait genetic variance,
- cross-trait genetic covariance,
- sum of squared realized marker effects.

It also sums across categories to get block-level total genetic variance.

#### 5. Sample residual covariance `R_blk`

If residual variance estimation is enabled:

1. The script samples a block-specific covariance matrix from an inverse-Wishart distribution.
2. It applies a heuristic threshold based on `sum(ssq_blk_cat) / totalvarg_blk` to decide whether to keep the sampled diagonal element or replace it with `1.0`.
3. It reconstructs the off-diagonal covariance from the sampled correlation and the retained diagonals.

After all blocks are processed, the script averages `R_blk` across blocks and reuses that common mean for the next iteration.

#### 6. Update derived quantities after burn-in and thinning

When `do_thin` is true, the script updates running posterior means for:

- `meanAlpha`
- marker-effect covariance estimated only from pleiotropic markers (`meanA`)
- sampled beta covariance (`meanB`)
- annotation-specific genetic covariance (`meanG`)
- total genetic covariance (`meanGtotal`)
- residual covariance (`meanR`)

The code keeps both first and second moments for many of these quantities so it can write posterior standard deviations at the end.

#### 7. Update mixture probabilities `Pi`

For each category, the script samples a four-state Dirichlet posterior using the state counts accumulated over all SNPs in that iteration.

The sampled result becomes the new `Pi[cat]`. If thinning is active, posterior means and second moments for `Pi` are updated as well.

#### 8. Update covariance hyperparameters

The code computes `SSE_vec` from the current `betaArray` values and uses it in two places.

- If `estimate_Gscale = true` (which can only be true when `estimate_vara = true`), during the initial `estGscale_iter` iterations it estimates `scale_G_vec` empirically from running averages of `beta'beta`.
- If `estimate_vara == true`, it samples each `A_vec[cat]` from an inverse-Wishart posterior with scale `scale_G_vec[cat] + SSE_vec[cat]`.

The inverse of each updated covariance is cached in `Ainv_vec` for the next marker update pass.

#### 9. Save checkpoints and sampled outputs

After burn-in, the script periodically writes:

- MCMC samples of `Pi`
- MCMC samples of beta covariance matrices
- MCMC samples of genetic covariance matrices
- sampled `delta` matrices at `outFreq`

At the final save point it also writes restart files, including:

- `last_mcmc_betaArray*.rank*.txt`
- `last_sample_R_blk/`
- `beta_effect_var_matrices_last_sample/`
- `pi_last_sample/`
- `last_sample_delta/`

### Final outputs

At the end of the run, the non-MPI workflow writes posterior summaries such as:

- `estA*.txt`: annotation-specific mean marker-effect covariance estimated from pleiotropic realized effects
- `estB*.txt`: annotation-specific mean sampled covariance of latent beta effects
- `estG*.txt`: annotation-specific mean genetic covariance
- `estGtotal.txt`: total genetic covariance across all annotations
- `estR.txt`: residual covariance
- `estPi*.txt`: annotation-specific posterior mean mixture probabilities
- `meanAlpha*.rank*.txt`: posterior mean realized SNP effects
- `mcmc_Delta*.rank*.txt`: saved inclusion trajectories

### Important implementation details

- The workflow accepts `seed`, but it does not currently call `Random.seed!`. The argument is currently informational unless that is added later.
- The workflow creates `analysis_path` before writing outputs.
- `Pi` state ordering is stabilized by helper functions such as `pi_key`, `pi_key_order`, and `write_pi_dict` rather than raw dictionary iteration order.

## Practical focus

### Use `scripts/run_nonmpi.jl` when

- you want the smallest runnable example,
- all transformed data fit into one process,
- you are debugging model logic or file conventions,
- you want to use the provided example under `example/`,
- you plan to add thread-level parallelism inside the current single-process workflow.

## Optimization-relevant observations

These points matter if the next step is code or repository optimization.

- The dominant compute cost is the nested loop over blocks, markers, categories, and traits.
- The code allocates many temporary vectors and slices inside hot loops, especially around `alphaArray`, `betaArray`, dictionary lookups, and per-marker trait updates.
- The next performance refactor should target thread-safe parallelism inside the current non-MPI workflow rather than reviving a second distributed entrypoint.
- The current restart and output format is text-heavy, which is convenient for inspection but expensive for large runs.
- The model state is stored in large dense vectors of length `my_nsnp * nCategory` for each trait, so memory pressure will rise quickly as the number of annotations grows.

That makes the current documentation step useful groundwork for later refactoring: most optimization work should preserve the statistical update order while reducing allocation, duplication, and synchronization overhead.