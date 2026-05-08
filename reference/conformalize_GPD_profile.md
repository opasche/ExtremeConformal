# GPD-profile conformalization procedure with optional sample weights

GPD-profile conformalization procedure with optional sample weights

## Usage

``` r
conformalize_GPD_profile(
  s_cal,
  coverage_proba,
  threshold_lvls,
  alpha_profile,
  alpha_profile_qlvl,
  init_step_pos = 100,
  init_step_neg = 10,
  tol = 0.001,
  steps_beyond_conf = 5,
  max_steps = 1000,
  obs_weights = NULL,
  verbose = 1,
  .pre_sorted = FALSE
)
```

## Arguments

- s_cal:

  .

- coverage_proba:

  .

- threshold_lvls:

  .

- alpha_profile:

  .

- alpha_profile_qlvl:

  .

- init_step_pos:

  .

- init_step_neg:

  .

- tol:

  .

- steps_beyond_conf:

  .

- max_steps:

  .

- obs_weights:

  .

- verbose:

  .

- .pre_sorted:

  .

## Value

A named list containing the conformal correction `dQ`, a redundant
`dQ_thresholds`, and placeholders for the estimated scale parameters
`sigma` and shape parameters `xi` (currently set to NULL).
