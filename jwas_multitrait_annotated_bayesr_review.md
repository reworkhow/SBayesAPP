# Review of JWAS Multi-Trait Annotated BayesC and BayesR

## Scope and source version

This note has two purposes:

1. audit SBayesAPP's current `annotation_prior_model = :marker_probit_tree`
   implementation against the documented JWAS multi-trait annotated BayesC
   algorithm; and
2. explain the newer dense two-trait annotated BayesR implementation on the
   JWAS `master` branch, as groundwork for a later SBayesAPP implementation.

The JWAS source review is pinned to commit
[`038c6c8a497ec6e17e957657a5c13ae698377132`](https://github.com/reworkhow/JWAS.jl/tree/038c6c8a497ec6e17e957657a5c13ae698377132),
dated 2026-07-11. This matters because the current
[multi-trait annotated BayesC manual](https://reworkhow.github.io/JWAS.jl/latest/manual/multitrait_annotated_bayesc/)
documents the three-step BayesC tree, while the separate annotated BayesR
manual still says that multi-trait BayesR is unsupported. The source code,
tests, design record, and benchmark reports at the pinned commit do contain a
dense two-trait annotated BayesR implementation. For that feature, the source
code is therefore the authoritative description.

## Bottom line

The current SBayesAPP `:marker_probit_tree` implementation is consistent with
the documented dense two-trait annotated BayesC model, specifically its
trait-at-a-time Gibbs sampler (JWAS sampler I). I found no mathematical
discrepancy in the four-state prior, latent-probit annotation update,
coefficient shrinkage prior, or marker-state update.

There are several implementation and interface differences, but they do not
change that conclusion:

- SBayesAPP stores the four states in the order `00, 11, 10, 01`; current JWAS
  stores them as `00, 10, 01, 11`. Both implementations map states explicitly,
  so this is a storage-order difference rather than a probability-model
  difference.
- SBayesAPP implements the trait-at-a-time sampler-I update. It does not expose
  JWAS's optional joint four-state sampler II.
- SBayesAPP operates on blockwise transformed summary-statistic inputs, whereas
  JWAS operates on individual-level genotypes.
- SBayesAPP constructs the starting four-state prior from single-trait results;
  JWAS accepts a joint starting `Pi` through its user API.
- SBayesAPP treats all supplied annotation columns as prior features in tree
  mode, forces `estimate_pi=true`, and currently disables continuation for this
  mode.

The newer JWAS feature is best understood as a separate annotated
multi-trait BayesR model, not merely as an extra mixture nested inside the
existing BayesC `11` state. It retains the BayesC `00/10/01/11` tree as its
first layer, then assigns a small, medium, or large BayesR magnitude class to
each active trait. This produces 16 joint states.

## Audit of the existing SBayesAPP tree

### Four joint inclusion states

For marker `j`, let

```text
d_j = (d_j1, d_j2), with d_jt in {0, 1}.
```

The possible states are:

- `00`: inactive for both traits
- `10`: active only for trait 1
- `01`: active only for trait 2
- `11`: active for both traits

SBayesAPP represents the corresponding latent and realized effects as

```text
beta_j  ~ N_2(0, A)
alpha_j = diag(d_j) * beta_j.
```

Thus an inactive trait still has a sampled latent `beta`, but its realized
effect `alpha` is zero. This is the same representation used by JWAS's
multi-trait BayesC sampler I.

### Three annotation probit steps

With annotation design row `a_j`, including an intercept, define

```text
p1_j = Phi(a_j' theta_1)
p2_j = Phi(a_j' theta_2)
p3_j = Phi(a_j' theta_3).
```

The steps are:

1. `p1_j = Pr(d_j != 00)`;
2. `p2_j = Pr(d_j = 11 | d_j != 00)`; and
3. `p3_j = Pr(d_j = 10 | d_j in {10, 01})`.

The marker-specific prior is reconstructed as

```text
Pr(00) = 1 - p1_j
Pr(11) = p1_j * p2_j
Pr(10) = p1_j * (1 - p2_j) * p3_j
Pr(01) = p1_j * (1 - p2_j) * (1 - p3_j).
```

This is implemented directly in
[`src/model/marker_probit_tree.jl`](src/model/marker_probit_tree.jl):

- `marker_probit_tree_step_indicators` builds the three binary responses and
  the correct conditional marker subsets;
- `sample_marker_probit_tree_step!` samples truncated-normal liabilities and
  regression coefficients; and
- `rebuild_marker_probit_tree_priors!` applies the equations above.

### Latent liabilities and coefficient priors

For each tree step `h`, SBayesAPP uses

```text
l_jh = a_j' theta_h + epsilon_jh,
epsilon_jh ~ N(0, 1),
z_jh = 1(l_jh > 0).
```

The latent probit residual variance is fixed at one. Liabilities are sampled
from the appropriate normal distribution truncated below or above zero.

The intercept is unpenalized. Each annotation slope has

```text
theta_hk | sigma_h^2 ~ N(0, sigma_h^2),  k > 1,
```

and SBayesAPP updates the step-specific shrinkage variance as

```text
sigma_h^2 = (sum_{k>1} theta_hk^2 + 2) / chi-square(p + 1),
```

where `p` is the total number of coefficients, including the intercept. This
matches the current JWAS implementation.

### Marker update

In
[`src/workflows/nonmpi.jl`](src/workflows/nonmpi.jl), SBayesAPP updates one
trait indicator at a time while holding the other trait's current state fixed.
For the current trait it evaluates the conditional Gaussian distribution and
the two possible joint states, reads their marker-specific probabilities from
`MarkerProbitTreeState.snp_pi`, samples the binary indicator, samples the
latent effect, and updates the working residual.

This is algebraically the same as JWAS's multi-trait BayesC sampler I. The
annotation model changes the joint prior term in this marker update; it does
not change the Gaussian likelihood calculation or create annotation-specific
effect-variance matrices.

### MCMC ordering

SBayesAPP uses the expected Gibbs order:

1. sample marker states and effects using the current marker-specific prior;
2. derive the three tree responses from the new marker states;
3. sample liabilities, annotation coefficients, and annotation shrinkage
   variances;
4. rebuild every marker's four-state prior for the next marker sweep; and
5. update the common latent marker-effect covariance matrix and the remaining
   variance components.

`Pi[1]` in SBayesAPP tree mode is only the genome-wide mean of the rows in
`snp_pi`, retained for summaries and output. Marker updates correctly use the
full marker-specific `snp_pi` matrix.

## What current JWAS adds: dense two-trait annotated BayesR

### It is a 16-state model

JWAS uses four BayesR class labels for each trait:

| JWAS class | Meaning | Variance multiplier `gamma` |
|---:|---|---:|
| 1 | zero | 0 |
| 2 | small | 0.01 |
| 3 | medium | 0.1 |
| 4 | large | 1.0 |

For marker `j`, the joint class state is

```text
c_j = (c_j1, c_j2), with c_jt in {1, 2, 3, 4}.
```

There are `4 x 4 = 16` states:

- one `00` state: `(1,1)`;
- three trait-1-only states: `(2:4,1)`;
- three trait-2-only states: `(1,2:4)`; and
- nine shared states: `(2:4,2:4)`.

Consequently, the BayesR magnitude mixture applies to every active trait, not
only to pleiotropic SNPs. A pleiotropic SNP can also have different magnitude
classes in the two traits, such as `(large, small)`.

### Layer 1: the original active-pattern tree

The first three probits are exactly the annotated BayesC tree:

```text
pattern00_j = 1 - p1_j
pattern11_j = p1_j * p2_j
pattern10_j = p1_j * (1 - p2_j) * p3_j
pattern01_j = p1_j * (1 - p2_j) * (1 - p3_j).
```

They determine whether the marker is inactive, trait-1-only, trait-2-only, or
shared.

### Layer 2: trait-specific BayesR magnitude trees

JWAS adds two probits per trait. For trait `t`, define

```text
q_t1,j = Pr(class_t is medium or large | trait t is active)
q_t2,j = Pr(class_t is large | class_t is medium or large).
```

Both probabilities are annotation-dependent probits. The conditional
probabilities of the three nonzero magnitude classes are

```text
m_t,small,j  = 1 - q_t1,j
m_t,medium,j = q_t1,j * (1 - q_t2,j)
m_t,large,j  = q_t1,j * q_t2,j.
```

Together with the three active-pattern steps, the annotation model therefore
has seven probit regressions:

1. zero versus any active pattern;
2. shared `11` versus a singleton pattern;
3. trait-1-only `10` versus trait-2-only `01`;
4. trait 1 medium-or-large versus small;
5. trait 1 large versus medium;
6. trait 2 medium-or-large versus small; and
7. trait 2 large versus medium.

All seven steps reuse the same standard-probit liability and coefficient
shrinkage machinery as the original tree.

### Reconstructing the full 16-state prior

For nonzero magnitude labels `k,l in {small, medium, large}`:

```text
Pr((0,0))   = pattern00_j
Pr((k,0))   = pattern10_j * m_1,k,j
Pr((0,l))   = pattern01_j * m_2,l,j
Pr((k,l))   = pattern11_j * m_1,k,j * m_2,l,j.
```

This factorization is an important modeling assumption: conditional on the
annotations and the shared/singleton pattern, the two traits' magnitude classes
factor into trait-specific probabilities. Cross-trait dependence is still
present through the probability of entering state `11` and through the marker
effect covariance matrix, but there is no additional annotation model for a
particular pair such as `(large, large)` beyond the product of the two
trait-specific magnitude probabilities.

The implementation is in JWAS
[`annotation_updates.jl`](https://github.com/reworkhow/JWAS.jl/blob/038c6c8a497ec6e17e957657a5c13ae698377132/src/1.JWAS/src/MCMC/annotation_updates.jl),
especially `bayesr_mt_step_indicators` and `rebuild_bayesr_mt_priors!`.

## Effect model and pleiotropy

JWAS samples a two-dimensional latent base effect

```text
beta_j ~ N_2(0, G).
```

For state `c_j`, define

```text
S_j = diag(sqrt(gamma[c_j1]), sqrt(gamma[c_j2])).
```

The realized effect used in the phenotype model is

```text
alpha_j = S_j * beta_j.
```

Therefore

```text
Var(alpha_j1 | c_j) = gamma[c_j1] * G11
Var(alpha_j2 | c_j) = gamma[c_j2] * G22
Cov(alpha_j1, alpha_j2 | c_j)
  = sqrt(gamma[c_j1] * gamma[c_j2]) * G12.
```

This is how the new model combines BayesR magnitudes with pleiotropy. `G12`
retains the cross-trait covariance of the latent effects, while the two class
multipliers scale the realized covariance for that marker.

## Marker-state sampler

Current JWAS does not sample all 16 states in one multinomial draw. It adapts
the existing trait-at-a-time multi-trait BayesC sampler I:

1. visit a marker;
2. update trait 1 while holding trait 2's current class fixed;
3. evaluate the four candidate classes for trait 1 using the corresponding
   entries in the marker's 16-state prior row;
4. sample the class and latent effect, then update the realized effect and
   residual; and
5. repeat for trait 2 using trait 1's newly sampled class.

For candidate class `r` of trait `k`, let

```text
s_r = sqrt(gamma[r])
s_o = sqrt(gamma[current class of the other trait])
q_j = x_j'x_j
Rinv = inv(R)
Ginv = inv(G).
```

The scalar conditional precision is

```text
C_r = Ginv[k,k] + q_j * s_r^2 * Rinv[k,k],
```

and the linear term is

```text
b_r = s_r * w_j' * Rinv[:,k]
      - (Ginv[k,o] + q_j * s_r * s_o * Rinv[k,o]) * beta_jo.
```

The conditional latent effect is

```text
beta_jk | r, rest ~ N(b_r / C_r, 1 / C_r),
```

and the class log weight is, up to a common constant,

```text
log_weight_r = log Pr(candidate joint state | annotations)
               - 0.5 * log(C_r)
               + 0.5 * b_r^2 / C_r.
```

Class 1, the zero class, still receives a latent `beta` draw. Its scale is zero,
so its realized `alpha` is zero. The implementation is in
[`MTBayesR.jl`](https://github.com/reworkhow/JWAS.jl/blob/038c6c8a497ec6e17e957657a5c13ae698377132/src/1.JWAS/src/markers/BayesianAlphabet/MTBayesR.jl).

## One JWAS MCMC iteration

For this path, the effective order is:

1. update fixed/location effects and working residuals;
2. sweep markers with the current marker-specific 16-state priors;
3. convert sampled class pairs into the seven binary annotation responses;
4. update the seven liability/probit regressions;
5. rebuild every marker's 16-state prior row;
6. replace the global `Pi` summary by the mean of those rows;
7. sample the global latent base-effect covariance `G`; and
8. sample the remaining variance components, including residual covariance.

The annotation probabilities are prior probabilities for the next marker
sweep. They are not themselves posterior inclusion probabilities. Posterior
trait PIPs and BayesR class frequencies are calculated from sampled class
assignments.

## Initialization and current support boundaries

At the pinned commit, the implemented multi-trait path has these properties:

- exactly two traits;
- dense genotype storage only;
- annotations are required;
- unconstrained full `G` covariance only;
- no random regression model;
- no streaming storage; and
- no fast-block multi-trait BayesR sampler.

The annotation state is initialized with zero coefficients, seven steps, 16
classes, and a uniform `1/16` prior row for every marker. Thus the first marker
sweep starts from a uniform joint state prior, after which the annotation
models learn marker-specific rows. This differs from single-trait annotated
BayesR, whose usual starting vector is `(0.95, 0.03, 0.015, 0.005)`.

The initialization and state ordering are defined in JWAS
[`annotation_setup.jl`](https://github.com/reworkhow/JWAS.jl/blob/038c6c8a497ec6e17e957657a5c13ae698377132/src/1.JWAS/src/markers/annotation_setup.jl).

## Covariance initialization and update

When JWAS converts the supplied genetic covariance to a starting marker-effect
covariance, it averages the appropriate BayesR scale products over the
marker-specific 16-state priors. Diagonal denominators use `gamma[class]` and
the off-diagonal denominator uses the product of square-root scales.

During MCMC, however, current JWAS updates `G` from

```text
sum_j beta_j * beta_j'
```

over all markers, including markers currently in the `00` state, and uses all
markers in the inverse-Wishart degrees of freedom. This is intentional in the
current implementation because the model samples latent `beta` for every
marker, but the JWAS implementation record explicitly lists an active-only
alternative as a future benchmark candidate. This choice should be treated as
a deliberate design decision when porting the method, not copied
unquestioningly.

See
[`variance_components.jl`](https://github.com/reworkhow/JWAS.jl/blob/038c6c8a497ec6e17e957657a5c13ae698377132/src/1.JWAS/src/variance_components.jl)
and the
[`multi-trait annotated BayesR implementation record`](https://github.com/reworkhow/JWAS.jl/blob/038c6c8a497ec6e17e957657a5c13ae698377132/docs/plans/2026-06-24-multitrait-annotated-bayesr-implementation.md).

## Interpretation of outputs

JWAS reports, for each trait and marker:

- posterior mean realized effect;
- posterior standard deviation;
- `Model_Frequency = Pr(class > 1 | data)`;
- separate posterior frequencies for classes 1 through 4; and
- `MediumLarge_Frequency = Pr(class >= 3 | data)`.

It also reports posterior summaries for all seven annotation coefficient
vectors. The current benchmark helper derives posterior `00/10/01/11`
probabilities from paired trait marker-effect samples; those collapsed pattern
probabilities do not retain the full 16-state magnitude information.

## Evidence from current JWAS benchmarks

JWAS includes a short smoke test and a longer production-quality benchmark for
this path. The benchmark authors' conclusion is cautious: the implementation
runs and is promising, but short chains produced unstable PIPs and too much
inactive-marker PIP in their simulated data. A 50,000-iteration run improved
average precision and cross-seed stability substantially, but the report did
not call the method fully production-ready for causal discovery.

This is relevant to a future SBayesAPP port: correctness tests should be
separated from mixing/calibration tests, and comparisons should use multiple
seeds and long enough chains.

Sources:

- [JWAS smoke-test report](https://github.com/reworkhow/JWAS.jl/blob/038c6c8a497ec6e17e957657a5c13ae698377132/benchmarks/reports/2026-06-24-multitrait-annotated-bayesr-simulated-annotations-report.md)
- [JWAS production-quality report](https://github.com/reworkhow/JWAS.jl/blob/038c6c8a497ec6e17e957657a5c13ae698377132/benchmarks/reports/2026-06-25-multitrait-annotated-bayesr-production-quality-report.md)

## Implications for a later SBayesAPP implementation

A clean extension should preserve the existing separation between the
annotation prior and the marker likelihood, but it will require more than
adding four columns to `MarkerProbitTreeState`.

The minimum coherent set of changes is likely to be:

1. add an explicit model or mixture-mode selector rather than silently changing
   existing `:marker_probit_tree` BayesC results;
2. replace binary per-trait indicators by class labels `1:4`;
3. expand marker-specific prior storage from four to 16 states;
4. expand the annotation state from three to seven probit steps;
5. change realized effects from `alpha = delta .* beta` to
   `alpha_t = sqrt(gamma[class_t]) * beta_t`;
6. replace the binary trait update by a four-class conditional update while
   retaining the full two-trait `R` and `G` cross terms;
7. initialize `G` with the expected BayesR scale products;
8. make an explicit, tested choice between JWAS's all-marker latent-`beta`
   update for `G` and an active-only alternative;
9. extend restart files, state-frequency outputs, and posterior summaries to
   retain all 16 states and class-specific frequencies; and
10. add unit tests for state indexing, row sums, seven response/active sets,
    zero-class behavior, conditional log weights, covariance scaling, and a
    short end-to-end run.

The existing SBayesAPP `:marker_probit_tree` path should remain a regression
baseline. Its current four-state behavior is internally coherent and matches
the documented annotated BayesC algorithm.
