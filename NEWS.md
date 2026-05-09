# KRLS 1.2-0

## New: formula interface

* `krls(y ~ x1 + x2, data = df)` is now supported alongside the
  matrix interface. The original `krls(X, y)` call continues to work
  unchanged. Drops the intercept column automatically; rejects NAs
  in either side with the same error messages the matrix interface
  uses.

## New: tidyverse-friendly extractors

* `tidy(fit)` returns a per-predictor data frame: `estimate` (the
  average marginal effect), `std.error`, `statistic`, `p.value`,
  `conf.low`/`conf.high`, plus pointwise quartiles
  (`q25`, `median`, `q75`) so users can see the heterogeneity that
  the AME averages over.
* `glance(fit)` is a one-row fit summary: `nobs`, `n_predictors`,
  `r.squared`, leave-one-out MSE, `lambda`, `sigma`, and effective
  degrees of freedom from the kernel ridge.
* `augment(fit, data = ...)` joins `.fitted`, `.resid`, and one
  `.dy_d_<predictor>` column per predictor (the pointwise
  derivatives) back to the data.

  Methods are registered against the `generics` package generics, so
  `library(broom)` makes them discoverable.

## New: `ggplot2` autoplot

* `autoplot(fit)` returns a faceted ggplot of the pointwise
  marginal-effect distribution, one panel per predictor, with the
  AME overlaid in blue. Discoverable via `library(ggplot2)`;
  `ggplot2` is in `Suggests:`.

## New: vignette

* `vignette("krls-quickstart", package = "KRLS")` walks through a
  small simulated example showing the AME-vs-pointwise heterogeneity
  story end-to-end.

# KRLS 1.1-0

## Bug fixes

* `lambdasearch()` now returns a plain scalar instead of a 1x1 matrix
  (it previously used `ifelse()` which inherits the shape of its test
  argument). This eliminates the R 4.4+ deprecation warning
  "Recycling array of length 1 in vector-array arithmetic is
  deprecated" that fired during `Eigenobject$values + lambda` in both
  `solveforc()` and `krls()`'s vcov calculation. The numerical results
  are identical to 1.0-0 within rounding.

* `plot.krls()` no longer dispatches `UseMethod("summary")` on a wrong-
  class input — that was a copy-paste mistake. All three S3 methods
  (`predict`, `summary`, `plot`) now `stop()` cleanly when the input is
  not a `krls` object, with a clear message. The class check itself is
  now `inherits(x, "krls")` rather than `class(x) != "krls"` (the old
  pattern misbehaves when the object has multiple classes).

## Internal cleanups (no user-visible change)

* `krls()`: removed dead code — pre-allocation of `Eigenobject$values`
  and `$vectors` immediately overwritten by `eigen()`, and several
  large commented-out alternative implementations.
* `solveforc()`: trimmed dead alternatives, modernized formatting.
* `predict.krls()`: removed dead commented `if(is.vector(newdata))`
  branch; fixed typo "standart errors" -> "standard errors".
* Redundant `sd(y) == 0` check in `krls()` removed (`var(y) == 0`
  already covers it).

## DESCRIPTION

* Switched to `Authors@R` (Jens Hainmueller as cre+aut, Chad Hazlett
  as aut) per current CRAN style.
* `Imports: grDevices, graphics, stats` added explicitly (these were
  already used via the NAMESPACE; making the dependency declaration
  explicit).
* `Suggests: testthat (>= 3.0.0)` added; `Config/testthat/edition: 3`
  added.
* `URL` field gains the GitHub repo and corrects the Stanford URL.
  `BugReports` field added.
* `Encoding: UTF-8` declared explicitly.
* Description text adds a DOI to the 2014 Political Analysis paper.
