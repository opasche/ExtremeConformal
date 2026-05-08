# GPD-max conformalization procedure with optional sample weights

GPD-max conformalization procedure with optional sample weights

## Usage

``` r
conformalize_GPD_max(
  s_cal,
  coverage_proba,
  threshold_lvls = NULL,
  min_obs_GPD = 10,
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

- min_obs_GPD:

  .

- obs_weights:

  .

- verbose:

  .

- .pre_sorted:

  .

## Value

A named list containing the conformal correction `dQ`, the vector of
`dQ_thresholds`, the estimated scale parameters `sigma`, the estimated
shape parameters `xi`, and the vector of effectively used thresholds
`threshold_lvls`.
