# cran-comments.md

## Submission notes for KRLS 1.1-0

This is a maintenance and modernization release relative to the
previously published KRLS 1.0-0 (2017-07-08).

### What's new

* **R 4.4+ deprecation fixed.** `lambdasearch()` previously returned a
  1x1 matrix (because `ifelse()` inherits the shape of its test
  argument). That matrix propagated into `Eigenobject$values + lambda`
  inside `solveforc()` and the vcov calculation in `krls()`,
  triggering "Recycling array of length 1 in vector-array arithmetic
  is deprecated" twice on every fit under recent R. The function now
  returns a plain scalar via `if (S1 < S2) X1 else X2`. Numerical
  results are unchanged.

* **`plot.krls()` and `predict.krls()` cleaner errors on bad input.**
  Both used `class(x) != "krls"` and (in `plot.krls`) dispatched
  `UseMethod("summary")` — a copy-paste error. All three S3 methods
  now use `inherits(x, "krls")` and `stop()` with a clear message.

* **DESCRIPTION modernized.** Switched to `Authors@R`, declared
  `Imports`, added `BugReports` and the GitHub `URL`, added `testthat`
  to `Suggests`, declared `Encoding: UTF-8`, added the 2014 PA paper
  DOI to the description text.

* **Several typo fixes** ("standart errors" -> "standard errors", etc.).

* **New tests/testthat/ suite** (34 assertions, 4 files) including a
  regression test for the R 4.4+ deprecation.

* **GitHub Actions CI workflow** added.

All changes preserve the byte-for-byte numerical results of 1.0-0 on
well-conditioned problems; verified by an internal regression test
(`dev/02_regression_check.R`) against the frozen 1.0-0 source.

### Test environments

* macOS Tahoe 26.3.1 (local), R 4.5.3
* Planned via win-builder R-devel and r-hub macOS-release on
  submission.

### R CMD check results

`R CMD check --as-cran` is clean locally — Status: 2 NOTEs, both
purely environmental and not present on CRAN's farm:

* "unable to verify current time" (network reach to a time server)
* "Skipping checking math rendering: package 'V8' unavailable"

Zero substantive NOTEs.

### Reverse dependencies

`KRLS` has 1 direct reverse dependency on CRAN: `InfluenceBorrowing`.
We ran `revdepcheck::revdep_check()` comparing 1.1-0 against the CRAN
baseline 1.0-0; no new problems were introduced.

### What we kept stable

* All exported function names and signatures: `krls`, `predict.krls`,
  `summary.krls`, `plot.krls`, `gausskernel`, `looloss`, `solveforc`.
* All `krls()` return-list field names. The function is fully
  backwards compatible.
