
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ExtremeConformal

<!-- badges: start -->

<!-- badges: end -->

This new extreme conformal prediction framework provides informative
prediction intervals at the high-confidence levels for which classical
conformal methods fail. In applications with potentially high-impact
events, a very high level of confidence is often required for
predictions. If that level is too large relative to the amount of data
used for calibration, classical conformal methods provide infinitely
wide, thus, uninformative prediction intervals. Our extreme conformal
procedure bridges extreme value statistics and conformal prediction to
provide reliable and informative prediction intervals with
high-confidence coverage, which can be constructed using any black-box
extreme quantile regression method. The methodology was introduced in
[Pasche, Lam, and Engelke
(2025)](https://doi.org/10.48550/arXiv.2505.08578)..

## Installation

To install the development version of ExtremeConformal from R, run

``` r
# install.packages("devtools")
devtools::install_github("opasche/ExtremeConformal")
```

------------------------------------------------------------------------

Package created by Olivier C. PASCHE  
Research Institute for Statistics and Information Science,  
University of Geneva (CH), 2025.  
Supported by the Swiss National Science Foundation.
