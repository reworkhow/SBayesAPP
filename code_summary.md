# SBayesAPP Code Summary

## What the code does

SBayesAPP implements a bivariate BayesC-style sampler that uses annotation-defined SNP groups to estimate:

- per-annotation mixture proportions for the four inclusion states across two traits,
- per-annotation SNP effect covariance matrices,
- per-annotation genetic covariance matrices,
- total genetic covariance matrices across the genome,
- residual covariance matrices, when that component is enabled.

The repository currently exposes two Julia entrypoints:

- `src/app_nonMPI.jl`: single-process version that assumes all LD blocks are stored in one input directory.
- `src/app_MPI.jl`: MPI-distributed version that uses the same core sampler, but splits LD blocks across ranks and synchronizes global parameters through MPI reductions and broadcasts.

The MPI code is not a different model. It is the same model with distributed data loading and distributed updates for global quantities.

## Shared model objects

Both scripts organize the sampler around the same parameterization.

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

The code stores those probabilities in dictionaries keyed by the corresponding length-2 vectors.

## Input contract

### Shared data inputs

Both scripts expect transformed block-wise inputs that have already been prepared upstream.

- `TransformedX_dict.jld2`: dictionary from block id to transformed genotype basis `X`.
- `TransformedY_dict.jld2`: dictionary from block id to transformed response vectors for trait 1 and trait 2.
- `blkSNPsIndex_dict.jld2`: per-block SNP indices.
- `blkIDs.txt`: block ids to analyze.
- `nGWAS_dict.jld2`: per-block sample size information.
- annotation matrix dictionary: a per-block annotation membership matrix saved as JLD2.
- annotation table: used to determine the annotation order and the number of SNPs in each annotation.

The current example input directory under `example/SBayesAPP_input_first10blks/` matches the non-MPI loader in `src/app_nonMPI.jl`.

### Single-trait initialization inputs

Both scripts set up several priors from a single-trait result directory `ST_path`, specifically:

- `Trait1/mean_pi.txt`
- `Trait2/mean_pi.txt`
- `Trait1/mean_varg_total.txt`
- `Trait2/mean_varg_total.txt`

These files are used to build:

- the initial mixture prior `startPi`,
- the initial scale of annotation-specific marker effect covariance matrix priors `Gprior_vec`.

### Annotation handling

Both scripts currently hard-code:

- `nCon = 0`
- `annotationType = repeat(["category"], nCat)`

That means the current checked-in code is operating in categorical-annotation mode only, even though helper logic exists for continuous annotations.

## `app_nonMPI.jl`

### Command-line interface

`src/app_nonMPI.jl` reads 13 positional arguments:

1. `data_path`
2. `analysis_path`
3. `nIter`
4. `seed`
5. `nrank`
6. `annot_file`
7. `annot_dict`
8. `outFreq`
9. `starting_value_dir`
10. `secondary_starting_value_dir`
11. `ST_path`
12. `thin`
13. `is_continue` (optional, default `false`)

In practice, `nrank` is expected to be `1` in this file.

### Setup phase

The single-process script does the following before entering MCMC:

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

The heart of the program is `runSBayesAPP`, whose iteration loop starts at `for iter = 1:nIter`.

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

At the end of the run, `app_nonMPI.jl` writes posterior summaries such as:

- `estA*.txt`: annotation-specific mean marker-effect covariance estimated from pleiotropic realized effects
- `estB*.txt`: annotation-specific mean sampled covariance of latent beta effects
- `estG*.txt`: annotation-specific mean genetic covariance
- `estGtotal.txt`: total genetic covariance across all annotations
- `estR.txt`: residual covariance
- `estPi*.txt`: annotation-specific posterior mean mixture probabilities
- `meanAlpha*.rank*.txt`: posterior mean realized SNP effects
- `mcmc_Delta*.rank*.txt`: saved inclusion trajectories

### Important implementation details

- The script parses `seed`, but unlike the MPI version it does not call `Random.seed!`. The argument is currently informational unless that is added later.
- The script assumes `analysis_path` already exists. The example shell script handles that with `mkdir -p`.
- Dictionary iteration order is relied on when matching the four `Pi` states to output rows. In modern Julia this is insertion-ordered, but it is still a fragile representation for a model parameter with fixed semantics.

## `app_MPI.jl`

### Command-line interface

`src/app_MPI.jl` extends the non-MPI argument list. It expects:

1. `data_path`
2. `analysis_path`
3. `nIter`
4. `seed`
5. `nrank`
6. `annot_file`
7. `annot_dict`
8. `outFreq`
9. `starting_value_dir`
10. `secondary_starting_value_dir`
11. `ST_path`
12. `thin`
13. `N1`
14. `N2`
15. `estimate_pi`
16. `fixed_hyperparameters` (optional)
17. `is_continue` (optional)

The extra arguments support:

- residual scaling based on average sample sizes `N1` and `N2`,
- turning `Pi` estimation on or off,
- fixing hyperparameters,
- continuing chains.

### Setup phase

The MPI script follows the same general preparation steps as the non-MPI version, with these additions.

#### MPI bootstrap

- `MPI.Init()` and `MPI.Finalize()` wrap the run.
- Each rank identifies `my_rank` and `cluster_size`.
- Random seeds are rank-specific via `Random.seed!(seed + my_rank)`.

#### Rank-local data loading

Each rank reads only its own shard from a directory structure of the form:

- `nrank<cluster_size>.eigen/bhatXsj/995Eigen/rank<my_rank>.*`

This is the key data-layout difference from the non-MPI script.

#### Optional execution modes

The MPI script also supports a fixed-hyperparameter mode.

- `fixed_hyperparameters`: disables updates for `Pi`, `A_vec`, and optionally `R`, and instead reads fixed values from prior outputs.

### How the sampler differs from `app_nonMPI.jl`

The local within-rank marker updates are almost the same as in `app_nonMPI.jl`. The main difference is what gets synchronized across ranks.

#### Quantities reduced to rank 0

At different stages of each iteration, the MPI script gathers and reduces to rank 0:

- `Atrue_vec` and pleiotropic marker counts `nQTL`
- annotation-state counts `nLoci_array_vec`
- annotation-specific genetic covariance `G_vec`
- total genetic covariance `G_total`
- beta sum-of-squares `SSE_vec`
- residual covariance totals `R_blk_sum`

These reductions let rank 0 update the global hyperparameters using statistics from the full data set rather than only its local shard.

#### Quantities broadcast from rank 0

After updating global parameters on rank 0, the script broadcasts:

- `Pi`
- `A_vec`

That keeps all ranks synchronized before the next marker-update pass.

#### Residual covariance handling

The residual covariance logic is slightly different from the non-MPI version.

- The prior starts at `diag(1/N1, 1/N2)`.
- Rank 0 averages residual covariance across all blocks using the total block count `nBlocks`.
- After reducing `R_blk_sum` and averaging it into `R_blkmean` on rank 0, the code currently does not broadcast `R_blkmean` back to all ranks. The relevant broadcast lines are commented out. As written, each rank keeps using its local `R_blk` values after the global summary is computed.

That point is an implementation detail worth remembering when comparing MPI and non-MPI behavior.

### Final outputs

Rank 0 writes global summaries analogous to the non-MPI outputs:

- `estA*.txt`, `estB*.txt`, `estG*.txt`, `estGtotal.txt`, `estR.txt`, `estPi*.txt`
- `mcmcGcov_*.txt`, `mcmcGcor_*.txt`, `mcmcAtruecor_c.txt`

Each rank also writes its own local posterior mean effects:

- `meanAlpha1.rank<rank>.txt`
- `meanAlpha2.rank<rank>.txt`

If delta saving is enabled, each rank also writes:

- `mcmc_Delta1.rank<rank>.txt`
- `mcmc_Delta2.rank<rank>.txt`

## Practical difference between the two entrypoints

### Use `app_nonMPI.jl` when

- you want the smallest runnable example,
- all transformed data fit into one process,
- you are debugging model logic or file conventions,
- you want to use the provided example under `example/`.

### Use `app_MPI.jl` when

- your transformed input has already been partitioned by rank,
- you want to run many LD blocks in parallel,
- you need a distributed version or fixed-hyperparameter workflows,
- you are running on a cluster with MPI available.

## Optimization-relevant observations

These points matter if the next step is code or repository optimization.

- The dominant compute cost is the nested loop over blocks, markers, categories, and traits.
- The code allocates many temporary vectors and slices inside hot loops, especially around `alphaArray`, `betaArray`, dictionary lookups, and per-marker trait updates.
- Both scripts duplicate a large amount of logic, so maintainability improvements should prioritize factoring out shared sampler code.
- The current restart and output format is text-heavy, which is convenient for inspection but expensive for large runs.
- The model state is stored in large dense vectors of length `my_nsnp * nCategory` for each trait, so memory pressure will rise quickly as the number of annotations grows.

That makes the current documentation step useful groundwork for later refactoring: most optimization work should preserve the statistical update order while reducing allocation, duplication, and synchronization overhead.