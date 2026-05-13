
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ExtremeConformal

<!-- badges: start -->

[![R-CMD-check](https://github.com/opasche/ExtremeConformal/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/opasche/ExtremeConformal/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This new extreme conformal prediction framework provides informative
prediction intervals at the high-confidence levels for which classical
conformal methods fail. Conformal prediction is a popular method to
construct prediction intervals with marginal coverage guarantees from
black-box machine learning models. In applications with potentially
high-impact events, a very high level of confidence is often required
for predictions. If that level is too large relative to the amount of
data used for calibration, classical conformal methods provide
infinitely wide, thus, uninformative prediction intervals. Our extreme
conformal procedure bridges extreme value statistics and conformal
prediction to provide reliable and informative prediction intervals with
high-confidence coverage, which can be constructed using any black-box
extreme quantile regression method. A weighted version of the approach
can account for nonstationary data. The methodology was introduced in
[Pasche, Lam, and Engelke
(2026)](https://doi.org/10.1007/s10687-026-00536-9).

## Installation

To install the development version of `ExtremeConformal`, in an R
session, run

``` r
# install.packages("devtools")
devtools::install_github("opasche/ExtremeConformal")
```

## References

Pasche, O. C., Lam, H., and Engelke, S. (2026). “Extreme Conformal
Prediction: Reliable Intervals for High-Impact Events.” *Extremes*.
[doi:10.1007/s10687-026-00536-9](https://doi.org/10.1007/s10687-026-00536-9).

------------------------------------------------------------------------

Package created by Olivier C. PASCHE  
Research Institute for Statistics and Information Science,  
University of Geneva (CH), 2025.  
Supported by the Swiss National Science Foundation.
