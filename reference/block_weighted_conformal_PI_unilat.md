# Single-sided block-weighted conformal prediction interval from conformalizer

Single-sided block-weighted conformal prediction interval from
conformalizer

## Usage

``` r
block_weighted_conformal_PI_unilat(
  Q_pred,
  block_ids,
  dQ_blocks,
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

- block_ids:

  Vector of block IDs for each test point, of the same length as
  `Q_pred`. The block IDs should be integers between 1 and the number of
  blocks, and correspond to the order of the `dQ_blocks` vector.

- dQ_blocks:

  Vector of block-specific conformal corrections to be added to the
  extreme quantile regression predictions, of the same length as the
  number of blocks.

- return_format:

  Format of the returned prediction interval. Either 'upper_limit' for a
  numerical vector of upper PI limits (default), 'interval' for a data
  frame with columns of lower and upper PI limits, or 'text' for a
  single string description of the PI (only available for single
  predictions).

- ymin:

  Lower endpoint for the response distribution (if known). Can be a
  single value (marginal lower endpoint), or a vector of the same length
  as `Q_pred` (conditional lower endpoint), or a vector of the same
  length as `Q_pred` (block-conditional lower endpoint). Default is
  `-Inf`.

- coverage_proba, coverage_alpha:

  (Optional) Marginal coverage probability (or level alpha) for the
  conformal prediction interval. Only one of `coverage_proba` or
  `coverage_alpha` must be provided, as
  `coverage_alpha = 1 - coverage_proba`. Only used for certain
  `return_format` options.

## Value

Depending on the `return_format` argument, either a numerical vector of
upper prediction interval (PI) limits, a data frame with lower and upper
PI limits, and the block IDs, as columns (and optionally coverage
probability and alpha), or (only available for compatibility) a single
string description of the PI (for single predictions only).
