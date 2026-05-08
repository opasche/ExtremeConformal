# Classical split conformalization procedure with optional sample weights

Classical split conformalization procedure with optional sample weights

## Usage

``` r
conformalize_Empirical(
  s_cal,
  coverage_proba,
  obs_weights = NULL,
  test_weight = NULL,
  verbose = 1,
  .pre_sorted = FALSE
)
```

## Arguments

- s_cal:

  .

- coverage_proba:

  .

- obs_weights:

  .

- test_weight:

  .

- verbose:

  .

- .pre_sorted:

  .

## Value

A named list containing the conformal correction `dQ` and a redundant
`dQ_thresholds`.
