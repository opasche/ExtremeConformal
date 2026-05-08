# Single-sided extreme conformal prediction

The main function of this package, performing single-sided extreme
conformal prediction from extreme quantile regression predictions. This
function computes the conformal correction score-quantile
\\\hat{q}^{e}\_{\alpha}\\ (or \\\hat{q}\_{\alpha}\\) to be added to the
extreme quantile regression prediction(s) to obtain the desired PIs. It
can optionally perform weighted conformal prediction for nonexchangeable
data (e.g. under distribution shift or drift), by providing
calibration-sample weights.

## Usage

``` r
conformalize_EQR_unilat(
  EQR_pred,
  y_calibr,
  coverage_proba = NULL,
  coverage_alpha = NULL,
  method = c("GPD_safeprofile", "GPD_boot", "GPD_profile", "Empirical",
    "GPD_profile_naive", "GPD_delta", "GPD_max", "GPD_simple", "Hill"),
  threshold_lvls = NULL,
  alpha_correction = c("Sidak", "Bonferroni"),
  correction_prop = 0.5,
  R = 1004,
  min_obs_GPD = 10,
  profile_init_step_pos = 100,
  profile_init_step_neg = 10,
  profile_tol = 0.001,
  profile_steps_beyond_conf = 5,
  profile_max_steps = 1000,
  alpha_profile_naive = 0.01,
  obs_weights = NULL,
  test_weight = NULL,
  return_q_list = FALSE,
  verbose = 1,
  .pre_sorted_scores = NULL
)
```

## Arguments

- EQR_pred:

  Vector of extreme quantile predictions for the calibration data,
  ideally at level `coverage_proba` (same length as `y_calibr`).

- y_calibr:

  Vector of observed response values of the calibration data (same
  length as `EQR_pred`).

- coverage_proba, coverage_alpha:

  Marginal coverage probability (or level alpha) for the conformal
  prediction interval. Only one of `coverage_proba` or `coverage_alpha`
  must be provided, as `coverage_alpha = 1 - coverage_proba`.

- method:

  Method for constructing the conformal prediction interval. Defaults to
  `'GPD_safeprofile'`. See the Details section for more information.

- threshold_lvls:

  Threshold probability level for the GPD-based methods. The `'GPD_max'`
  model instead requires a vector of threshold levels. Defaults to 0.95
  (or to a sequence from 0.8 to 0.99 with the `'GPD_max'` method).

- alpha_correction:

  Confidence correction method for the CI-based extreme conformal
  prediction intervals used to choose \\\alpha_1\\ and \\\alpha_2\\.
  Defaults to `'Sidak'`.

- correction_prop:

  Relative (pseudo-)proportion between the values of \\\alpha_1\\ and
  \\\alpha_2\\, for the CI-based extreme conformal prediction intervals.
  Defaults to 0.5 (equal values of \\\alpha_1\\ and \\\alpha_2\\).
  Larger values increase \\\alpha_1\\ and decrease \\\alpha_2\\.

- R:

  Number of bootstrap replicates for the `''GPD_boot'` method.

- min_obs_GPD:

  Minimum number of observations above the threshold for the GPD-based
  methods.

- profile_init_step_pos, profile_init_step_neg:

  Initial binary-search step size, in the positive and negative
  directions, for the profile-likelihood methods. See
  [`ExtremeCI::GPD_profile_CI()`](https://opasche.github.io/ExtremeCI/reference/GPD_profile_CI.html)
  for more details.

- profile_tol:

  Tolerance for the profile-likelihood search. See
  [`ExtremeCI::GPD_profile_CI()`](https://opasche.github.io/ExtremeCI/reference/GPD_profile_CI.html)
  for more details.

- profile_steps_beyond_conf:

  Number of safety steps beyond the confidence line for the
  profile-likelihood search. See
  [`ExtremeCI::GPD_profile_CI()`](https://opasche.github.io/ExtremeCI/reference/GPD_profile_CI.html)
  for more details.

- profile_max_steps:

  Maximum number of initial search steps for the profile-likelihood CIs.
  See
  [`ExtremeCI::GPD_profile_CI()`](https://opasche.github.io/ExtremeCI/reference/GPD_profile_CI.html)
  for more details.

- alpha_profile_naive:

  CI confidence level (i.e., the equivalent of \\\alpha_2\\) for the
  `'GPD_profile_naive'` method.

- obs_weights:

  Optional vector of sample weights (same length as `y_calibr`), to
  perform weighted conformal prediction for nonexchangeable data (e.g.
  under distribution shift or drift).

- test_weight:

  Weight of the test point for which the weighted conformal prediction
  is performed. Only necessary for the `'Empirical'` method if
  `obs_weights` is provided (defaults to `max(obs_weights)`).

- return_q_list:

  Boolean indicating whether to return a debug quantile list.

- verbose:

  Verbose level (0 for no messages, 1 for warnings, 2 for warnings
  duplicated as [`cat()`](https://rdrr.io/r/base/cat.html) prints).

- .pre_sorted_scores:

  (For development only.) If already computed: Vector of pre-sorted
  calibration nonconformity scores.

## Value

A named list containing the following elements.

- dQ:

  The conformal correction \\\hat{q}^{e}\_{\alpha}\\ or
  \\\hat{q}\_{\alpha}\\ to be added to the extreme quantile
  prediction(s) to obtain the conformal PI endpoint.

- coverage_proba:

  The marginal coverage probability for the conformal PI.

- coverage_alpha:

  The marginal coverage alpha level for the conformal PI.

- threshold_lvls:

  The threshold level(s) effectively used for the GPD-based methods.

- method:

  The `method` used to obtain the conformal PI.

- dQ_thresholds:

  A debugging vector of quantiles at the specified threshold levels.
  Only for GPD-based methods when `return_q_list==TRUE`.

- sigma:

  The GPD scale parameter estimate. Currently only for method
  'GPD_simple', when `return_q_list==TRUE`.

- xi:

  The GPD shape parameter estimate. Currently only for method
  'GPD_simple', when `return_q_list==TRUE`.

## Details

The `method` argument specifies the conformalization method used to
construct the (extreme) conformal prediction intervals (PIs). The
available method options are:

- 'GPD_safeprofile':

  Recommended choice for extreme conformal prediction. Tries the method
  'GPD_profile' first, and falls back to 'GPD_boot' if the former
  suffers from numerical instability.

- 'GPD_profile':

  GPD-based extreme conformalization using the profile-likelihood CI
  endpoint of the extreme score quantile. It captures the asymetric
  uncertainty of the score-quantile best, yielding the most reliable
  coverage. It might sometimes overcover or suffer from numerical
  convergence issues.

- 'GPD_boot':

  GPD-based extreme conformalization using the nonparametric bootstrap
  percentile CI endpoint of the extreme score quantile.

- 'GPD_delta':

  GPD-based extreme conformalization using the Delta method CI endpoint
  of the extreme score quantile.

- 'Empirical':

  The classical (non-extreme) conformalized quantile regression method,
  relying on the empirical quantile of the nonconformity scores. Yields
  infinitely wide PIs if the `coverage_proba` is larger than
  `1-1/(length(y_calibr)+1`).

- 'GPD_profile_naive':

  A naive version of `'GPD_profile'`, without the `alpha_correction` for
  multiple testing. Is likely to undercover, only use for comparison.

- 'GPD_max':

  A naive GPD-based approach repeating the `'GPD_simple'` method for a
  range of threshold levels, keeping the most conservative results. Is
  likely to undercover, only use for comparison.

- 'GPD_simple':

  A naive GPD-based extreme conformalization using a simple extrapolated
  score quantile estimate instead of a CI endpoint. Is likely to
  undercover, only use for comparison.

- 'Hill':

  (Not implemented) Extreme conformalization based on the Hill estimator
  from extreme value analysis.

See Pasche et al. (2026), referenced below, for the technical details of
extreme conformal prediction.

## References

Pasche, O. C., Lam, H., and Engelke, S. (2026). "Extreme Conformal
Prediction: Reliable Intervals for High-Impact Events." *Extremes*.
[doi:10.1007/s10687-026-00536-9](https://doi.org/10.1007/s10687-026-00536-9)
.
