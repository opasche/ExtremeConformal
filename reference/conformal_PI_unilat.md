# Single-sided conformal prediction interval from conformalizer

Single-sided conformal prediction interval from conformalizer

## Usage

``` r
conformal_PI_unilat(
  Q_pred,
  unilat_conformalizer,
  return_format = c("upper_limit", "interval", "text"),
  ymin = c(-Inf),
  coverage_proba = NULL,
  coverage_alpha = NULL
)
```

## Arguments

- Q_pred:

  Vector of extreme quantile regression predictions for the test data
  (same length as the number of test points).

- unilat_conformalizer:

  Either a conformalizer object obtained from
  [`conformalize_EQR_unilat()`](https://opasche.github.io/ExtremeConformal/reference/conformalize_EQR_unilat.md),
  or, directly, the conformal correction to be added to the extreme
  quantile regression predictions. In the latter case, either a single
  conformal correction value or a vector of conformal corrections of the
  same length as `Q_pred` is expected.

- return_format:

  Format of the returned prediction interval. Either 'upper_limit' for a
  numerical vector of upper PI limits (default), 'interval' for a data
  frame with columns of lower and upper PI limits, or 'text' for a
  single string description of the PI (only available for single
  predictions).

- ymin:

  Lower endpoint for the response distribution (if known). Can be a
  single value (marginal lower endpoint) or a vector of the same length
  as `Q_pred` (conditional lower endpoint). Default is `-Inf`.

- coverage_proba, coverage_alpha:

  (Optional) Marginal coverage probability (or level alpha) for the
  conformal prediction interval. Only one of `coverage_proba` or
  `coverage_alpha` must be provided, as
  `coverage_alpha = 1 - coverage_proba`. Only used for certain
  `return_format` options.

## Value

Depending on the `return_format` argument, either a numerical vector of
upper prediction interval (PI) limits, a data frame with lower and upper
PI limits as columns (and optionally coverage probability and alpha), or
a single string description of the PI (for single predictions only).

## Examples

``` r
conformal_PI_unilat(Q_pred=c(10,12), unilat_conformalizer=0.3, return_format='upper_limit')
#> [1] 10.3 12.3
```
