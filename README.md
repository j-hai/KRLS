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

## Scaling to large n: Nyström approximation

For samples where the exact $n \times n$ kernel matrix is
uncomfortable (typically past $n \approx 5{,}000$), KRLS offers an
explicit low-rank alternative since version 1.4-0:

```r
fit <- krls(X = X, y = y, approx = "nystrom")
```

Same interface, same fit object, same `summary()` / `predict()` /
`tidy()` outputs — including standard errors for predictions and
average marginal effects. Time becomes $O(n m^2 + m^3)$ and memory
$O(n m)$ with $m = \min(500, \lceil\sqrt{n}\rceil)$ landmarks by
default. See `vignette("krls-nystrom-scaling")` for a timing
comparison, the landmark-reuse pattern via `get_landmarks()`, and
the LOO-vs-GCV trade-off for selecting $\lambda$.

## What's new since 1.1-0

* **`approx = "nystrom"` (1.4-0)** — explicit low-rank fit path
  with analytical SEs for predictions and AMEs. Random and k-means
  landmark selection.
* **AME-variance optimization (1.4-0)** — row-sum rewrite reduces
  per-predictor work from $O(n^3)$ to $O(n^2)$; default `krls()`
  is 1.2–3× faster at moderate $n$, bit-identical results.
* **`get_landmarks()` accessor and diagnostic-rich `summary()` /
  `glance()` (1.4-1)** — extract landmarks for sensitivity-check
  refits without the standardize-twice gotcha; surface
  approximation diagnostics.
* **`lambda_method = "gcv"` (1.5-0)** — closed-form GCV alternative
  to the historical leave-one-out criterion. Default stays
  `"loo"`.
* **Nyström scaling vignette (1.5-0)** — exact-vs-Nyström timing,
  landmark reuse, LOO-vs-GCV.
* **Bug fixes (1.4-1, 1.5-1, 1.5-2)** — `tol` propagation,
  eigtrunc single-eigenvector path, Nyström at $m = 1$, k-means
  with $m = n$, GCV print label, clear error when the cross-kernel
  underflows.

See [`NEWS.md`](NEWS.md) for the full change log.

## License

GPL (>= 2). See `LICENSE.note`.
