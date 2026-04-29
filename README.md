# KRLS

[![R-CMD-check](https://github.com/j-hai/KRLS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/j-hai/KRLS/actions/workflows/R-CMD-check.yaml)

R implementation of **Kernel-Based Regularized Least Squares (KRLS)** —
a machine-learning method for regression and classification that fits a
flexible function `y = f(x)` without assuming linearity or additivity.
KRLS solves a Tikhonov regularization problem with a squared loss using
Gaussian kernels as radial basis functions, yielding closed-form
estimates of fitted values, variances, and pointwise marginal effects
that can be inspected for heterogeneity and interactions.

The method is described in:

> Hainmueller, J., & Hazlett, C. (2014). "Kernel Regularized Least Squares: Reducing Misspecification Bias with a Flexible and Interpretable Machine Learning Approach." *Political Analysis*, 22(2), 143–168. [doi:10.1093/pan/mpt019](https://doi.org/10.1093/pan/mpt019)

## Installation

```r
# From CRAN (stable)
install.packages("KRLS")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("j-hai/KRLS")
```

## Quick start

```r
library(KRLS)

# Toy data with known linear truth
set.seed(1)
N  <- 200
x1 <- rnorm(N)
x2 <- rbinom(N, 1, 0.3)
y  <- x1 + 0.5 * x2 + rnorm(N, 0, 0.2)
X  <- cbind(x1, x2)

# Fit
fit <- krls(X = X, y = y)

# Inspect — marginal effects with standard errors and t-stats
summary(fit)

# Visualize — distribution of pointwise marginal effects, plus
# conditional expectation curves for each predictor
plot(fit)

# Predict on new data
Xnew <- cbind(rnorm(20), rbinom(20, 1, 0.3))
pr   <- predict(fit, newdata = Xnew, se.fit = TRUE)
pr$fit
pr$se.fit
```

## What's new in 1.1-0

* Fixed an R 4.4+ deprecation that fired twice on every fit
  (`Eigenobject$values + lambda` recycled a 1×1 matrix from
  `lambdasearch()`). Numerical results are unchanged.
* `plot.krls()` no longer dispatches the wrong S3 generic on bad
  input. All three S3 methods now `stop()` cleanly and use
  `inherits()` for the class check.
* DESCRIPTION modernized (`Authors@R`, `Imports`, `Suggests`,
  `BugReports`, GitHub `URL`, DOI in description).
* Several typos fixed in error messages.
* New `tests/testthat/` suite covering `krls()`, `predict()`,
  `summary()`, helpers, and the deprecation regression.
* GitHub Actions CI on macOS / Windows / Ubuntu × release + devel
  + oldrel-1.

See [`NEWS.md`](NEWS.md) for the full change log.

## License

GPL (>= 2). See `LICENSE.note`.
