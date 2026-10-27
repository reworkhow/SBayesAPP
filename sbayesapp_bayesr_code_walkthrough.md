# SBayesAPP Two-Trait Annotated BayesR Code Walkthrough

This guide follows the order in which the program actually runs. The important
distinction is that `src/workflows/nonmpi.jl` is the driver, while
`src/model/marker_probit_tree.jl` contains helper functions called by that
driver. Reading `marker_probit_tree.jl` from top to bottom does not reproduce
the runtime order.

The BayesR mode is selected with:

```text
annotation_prior_model = :marker_probit_tree_bayesr
```

## One-screen call tree

```text
scripts/run_nonmpi.jl
└── parse_nonmpi_args                         src/SBayesAPP.jl:140
    └── NonMPIConfig                          src/config/types.jl:84
└── run_nonmpi                               src/SBayesAPP.jl:197
    └── run_nonmpi_workflow                  src/workflows/nonmpi.jl:1012
        ├── build_nonmpi_run_context         src/workflows/nonmpi.jl:324
        │   ├── build_nonmpi_sampler_settings
        │   ├── load_annotation_metadata
        │   ├── build_start_pi
        │   └── defer BayesR G initialization until the SNP count is known
        └── run_nonmpi_sampler!              src/workflows/nonmpi.jl:389
            ├── load_nonmpi_block_data
            ├── build_bayesr_gprior_vec
            ├── initialize_nonmpi_parameter_state
            ├── prepare_marker_probit_tree_block_state!
            ├── initialize_marker_probit_tree_state
            ├── initialize beta, alpha, and delta
            └── for each MCMC iteration
                ├── run_nonmpi_block_range!  sample marker classes/effects
                ├── update_marker_probit_tree_priors!
                ├── update latent beta covariance
                ├── summarize realized genetic covariance
                └── save retained samples
            └── save_nonmpi_posterior_mean!
```

## Important state variables

| Variable | Meaning in BayesR mode |
|---|---|
| `deltaArray[t]` | BayesR class `1:4` for trait `t`, not a binary indicator |
| `betaArray[t]` | latent base effect for trait `t`; sampled even when class is 1 |
| `alphaArray[t]` | realized effect, `sqrt(gamma[class]) * beta` |
| `gamma` | `(0, 0.01, 0.1, 1.0)` for classes `1:4` |
| `marker_probit_tree_state.snp_pi` | marker-by-16 matrix of annotation-induced joint class priors |
| `A_vec[1]` | common covariance matrix of the latent two-trait `beta` effects |
| `Ainv_vec[1]` | inverse of `A_vec[1]`, used during marker sampling |
| `R_blk[b]` | residual covariance matrix for LD block `b` |
| `Pi[1]` | genome-average `00/11/10/01` summary; marker sampling uses `snp_pi`, not this summary |

## Phase 1: load the package and parse arguments

1. `scripts/run_nonmpi.jl:1-7` activates the project, loads `SBayesAPP`, parses
   command-line arguments, and calls `run_nonmpi`.
2. `src/SBayesAPP.jl:3-18` includes all source files. This is source-loading
   order, not sampler execution order.
3. `parse_nonmpi_args` at `src/SBayesAPP.jl:140-187` reads
   `--annotation_prior_model marker_probit_tree_bayesr` and constructs a
   `NonMPIConfig`.
4. `src/config/types.jl:5-9` recognizes the BayesR mode as valid.
5. `NonMPIConfig` at `src/config/types.jl:39-82` stores the run options.
   `MarkerProbitTreeState` at lines 68-82 defines storage for the annotation
   design, seven probit regressions, liabilities, and the 16-state prior.
6. `build_nonmpi_sampler_settings` at
   `src/config/nonmpi_settings.jl:1-31` applies tree-mode rules: `n_con` is
   ignored, `estimate_pi` is forced on, and continuation is disabled.

## Phase 2: build the run context and read input

`run_nonmpi_workflow` starts at `src/workflows/nonmpi.jl:1012`:

1. It sets the random seed and restricts BLAS to one thread when Julia block
   threading is enabled.
2. It calls `build_nonmpi_run_context` at
   `src/workflows/nonmpi.jl:324-340`.
3. `load_annotation_metadata` at `src/io/inputs.jl:59-77` reads the annotation
   table and records annotation names. In a marker-tree model, every annotation
   column is treated as a prior predictor.
4. `build_start_pi` in `src/model/priors.jl:5-24` reads the single-trait
   starting results. In BayesR mode this four-pattern value is retained for
   compatible summaries, but the marker-specific 16-state prior starts
   uniformly.
5. `run_nonmpi_sampler!` begins at `src/workflows/nonmpi.jl:389`.
6. `load_nonmpi_block_data` is called at line 416 and implemented at
   `src/io/inputs.jl:403-447`. It reads transformed `X`, transformed `Y`, SNP
   indices, block IDs, GWAS sample sizes, and block annotation matrices into
   `NonMPIBlockData` objects.

## Phase 3: initialize the BayesR sampler

Initialization occurs inside `run_nonmpi_sampler!` before the MCMC loop:

1. `src/workflows/nonmpi.jl:431-438` sets `nCategory = 1` and calls
   `build_bayesr_gprior_vec`.
2. `build_bayesr_gprior_vec` at `src/model/priors.jl:40-52` converts the
   single-trait genetic variances into a starting latent-`beta` covariance by
   dividing by the expected BayesR multiplier under the uniform start.
3. `initialize_nonmpi_parameter_state` is called at
   `src/workflows/nonmpi.jl:441-456` and implemented at
   `src/model/initial_state.jl:4-76`. It creates `A_vec`, `Ainv_vec`, the
   compatibility `Pi` summary, and priors for `R` and `A`.
4. `prepare_marker_probit_tree_block_state!` is called at
   `src/workflows/nonmpi.jl:468` and implemented at
   `src/model/block_setup.jl:112-120`. It creates one all-marker effect
   category and one global annotation design matrix with an intercept.
5. `initialize_marker_probit_tree_state` is called at
   `src/workflows/nonmpi.jl:469-473` and implemented at
   `src/model/marker_probit_tree.jl:71-100`. For BayesR it allocates seven
   probit steps, 16 prior columns, zero coefficients, and a uniform `1/16`
   starting prior for every marker.
6. `src/workflows/nonmpi.jl:490-500` allocates compact `betaArray`,
   `alphaArray`, and `deltaArray`. Every BayesR class is initialized to class 1,
   while both latent and realized effects start at zero.
7. `src/workflows/nonmpi.jl:518-547` allocates posterior summaries and output
   files. If `output_mcmc_delta=true`, space is reserved for retained class
   labels.

## Phase 4: start one MCMC iteration

The loop begins at `src/workflows/nonmpi.jl:619`.

The effective order for each iteration is:

```text
1. sample marker classes, latent beta, and realized alpha
2. update the seven annotation probit models from the sampled classes
3. rebuild each marker's 16-state prior for the next iteration
4. update the common latent-beta covariance matrix
5. summarize residual and realized genetic covariance
6. save the iteration if it is retained
```

The marker sweep therefore uses the `snp_pi` produced at the end of the
previous iteration. The first sweep uses the uniform `1/16` initialization.

## Phase 4A: sample marker classes and effects

At `src/workflows/nonmpi.jl:639-708`, the driver calls
`run_nonmpi_block_range!`, either directly or once per thread.

Inside `run_nonmpi_block_range!`:

1. `src/workflows/nonmpi.jl:112-145` visits each assigned LD block, calculates
   `Rinv`, and randomizes marker order.
2. Lines 147-174 visit one marker and load its current `beta`, `alpha`, class
   label, and partial-residual crossproduct.
3. Lines 176-201 enter the BayesR-specific branch.
4. Trait 1 is updated first while holding trait 2's current class and latent
   effect fixed. Trait 2 is then updated using trait 1's newly sampled values.
5. `bayesr_conditional_class_parameters!` at
   `src/model/marker_probit_tree.jl:348-379` calculates the conditional
   precision, mean, variance, and log weight for each candidate class `1:4`.
   The log prior comes from that marker's row in `snp_pi`.
6. `sample_log_categorical!` at
   `src/model/marker_probit_tree.jl:381-396` normalizes the four log weights
   safely and samples the class.
7. `src/workflows/nonmpi.jl:197-200` samples latent `beta`, constructs realized
   `alpha = sqrt(gamma[class]) * beta`, and updates the working residual.
8. Lines 247-254 write the new `beta`, class label, and `alpha` back into the
   compact arrays.
9. Lines 259-302 calculate block-level realized genetic-variance summaries
   using `alpha`.
10. Lines 304-320 optionally sample the block residual covariance `R_blk[b]`.

Class 1 has `sqrt(gamma[1]) = 0`, so its realized `alpha` is zero. Its latent
`beta` is nevertheless sampled and later contributes to the latent covariance
update.

## Phase 4B: update the annotation prior

After all blocks have completed their marker sweeps, execution reaches
`src/workflows/nonmpi.jl:749-785`:

1. Line 772 expands the compact class arrays into marker order.
2. Line 773 calls `update_marker_probit_tree_priors!` at
   `src/model/marker_probit_tree.jl:325-339`.
3. `marker_probit_tree_bayesr_step_indicators` at
   `src/model/marker_probit_tree.jl:130-161` converts the new class pairs into
   seven binary responses and seven appropriate active subsets.
4. For each of the seven steps, `sample_marker_probit_tree_step!` at
   `src/model/marker_probit_tree.jl:220-257` performs one probit update.
5. That step calls `sample_marker_probit_tree_liabilities!` at lines 163-181 to
   sample truncated-normal liabilities according to the binary response.
6. It calls `update_marker_probit_tree_coefficients!` at lines 183-210 to
   sample the intercept and annotation coefficients.
7. It calls `sample_marker_probit_tree_step_variance!` at lines 212-218 to
   update the annotation-slope shrinkage variance.
8. After all seven steps,
   `rebuild_marker_probit_tree_bayesr_priors!` at lines 271-297 converts the
   seven probit probabilities into all 16 joint class probabilities for every
   marker. These priors will be used by the next iteration's marker sweep.
9. `marker_probit_tree_summary_dict` at lines 299-314 collapses the 16-state
   prior into genome-average `00/11/10/01` values for legacy `Pi` output only.
10. For retained iterations,
    `record_marker_probit_tree_coefficient_moments!` at lines 316-323 updates
    running posterior moments of coefficients and 16-state prior summaries.

## Phase 4C: update latent effect covariance

`src/workflows/nonmpi.jl:787-841` updates the common covariance of latent
effects:

1. Lines 791-801 calculate `SSE_vec = beta' * beta` using all markers,
   including class-1 markers.
2. Lines 803-816 optionally adapt the inverse-Wishart scale during the initial
   `estGscale_iter` iterations.
3. Lines 818-834 sample `A_vec[1]` from its inverse-Wishart conditional.
4. Lines 836-841 calculate `Ainv_vec[1]` for the next marker sweep.

Do not confuse this latent covariance `A_vec` with the realized genetic
covariance calculated from `alpha` later in the iteration.

## Phase 4D: summarize and save

1. `src/workflows/nonmpi.jl:843-854` accumulates residual covariance summaries.
2. Lines 856-894 calculate realized genetic covariance from `X * alpha`.
3. Lines 896-947 save retained `Pi`, covariance, class, and restart values.
4. After the MCMC loop, lines 955-995 call
   `save_nonmpi_posterior_mean!`.
5. `save_nonmpi_posterior_mean!` at `src/io/outputs.jl:125-263` writes the
   posterior means, seven-step annotation coefficients, and
   `annotation_bayesr_joint_prior.txt`.
6. `src/workflows/nonmpi.jl:1000-1003` writes `mcmc_Delta1.txt` and
   `mcmc_Delta2.txt`. In BayesR mode these contain class labels `1:4`.

## Functions in `marker_probit_tree.jl`, in runtime-use order

The functions are defined in a convenient source organization, but their
BayesR runtime order is:

| Runtime stage | Function | Purpose |
|---|---|---|
| Initialization | `initialize_marker_probit_tree_state` | allocate seven-step and 16-state annotation state |
| Marker sweep | `bayesr_conditional_class_parameters!` | calculate four class conditionals for one trait |
| Marker sweep | `log_marker_probit_tree_state_prior` | retrieve the selected joint-state prior for one marker |
| Marker sweep | `marker_probit_tree_bayesr_state_index` | map `(class1,class2)` to one of 16 columns |
| Marker sweep | `marker_probit_tree_bayesr_scale` | return `sqrt(gamma[class])` |
| Marker sweep | `sample_log_categorical!` | sample one of the four candidate classes |
| Annotation update | `update_marker_probit_tree_priors!` | coordinate the entire seven-step update |
| Annotation update | `marker_probit_tree_bayesr_step_indicators` | derive the seven binary responses and subsets |
| Each probit step | `sample_marker_probit_tree_step!` | coordinate liability, coefficient, and variance updates |
| Each probit step | `sample_marker_probit_tree_liabilities!` | sample truncated-normal liabilities |
| Each probit step | `update_marker_probit_tree_coefficients!` | sample probit regression coefficients |
| Each probit step | `sample_marker_probit_tree_step_variance!` | sample slope shrinkage variance |
| Annotation update | `rebuild_marker_probit_tree_bayesr_priors!` | reconstruct all marker-by-16 prior probabilities |
| Summary | `marker_probit_tree_summary_dict` | collapse 16 states to the four activity patterns |
| Retained sample | `record_marker_probit_tree_coefficient_moments!` | update posterior moments |

## Recommended reading order

For a line-by-line code review, read only the indicated sections in this order:

1. `src/config/types.jl:5-9,39-82`
2. `scripts/run_nonmpi.jl:1-7`
3. `src/SBayesAPP.jl:140-197`
4. `src/workflows/nonmpi.jl:1012-1018`
5. `src/workflows/nonmpi.jl:324-340`
6. `src/io/inputs.jl:59-77,403-447`
7. `src/workflows/nonmpi.jl:389-619`
8. `src/workflows/nonmpi.jl:62-323`
9. `src/model/marker_probit_tree.jl:348-396`
10. `src/workflows/nonmpi.jl:749-841`
11. `src/model/marker_probit_tree.jl:130-161,163-257,271-339`
12. `src/workflows/nonmpi.jl:843-1009`
13. `src/io/outputs.jl:125-263`

The mathematical tests in `test/marker_probit_tree_bayesr.jl` provide a short
executable specification of the 16-state indexing, seven responses, prior row
sums, collapsed probabilities, and zero-class conditional behavior.
