# GPD-boot conformalization procedure with optional sample weights

GPD-boot conformalization procedure with optional sample weights

## Usage

``` r
conformalize_GPD_boot(
  s_cal,
  threshold_lvls,
  alpha_profile,
  alpha_profile_qlvl,
  R = 1004,
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

- R:

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
