# GPD-delta conformalization procedure with optional sample weights

GPD-delta conformalization procedure with optional sample weights

## Usage

``` r
conformalize_GPD_delta(
  s_cal,
  threshold_lvls,
  alpha_profile,
  alpha_profile_qlvl,
  obs_weights = NULL,
  verbose = 1,
  .pre_sorted = FALSE
)
```

## Arguments

- s_cal:

  .

- threshold_lvls:

  .

- alpha_profile:

  .

- alpha_profile_qlvl:

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
