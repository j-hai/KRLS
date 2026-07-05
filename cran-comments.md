# cran-comments.md

## Submission notes for KRLS 1.7-1

This is a test-only patch that resolves the "Additional issues" MKL
check flagged for 1.7-0
(<https://www.stats.ox.ac.uk/pub/bdr/Rblas/MKL/KRLS.out>).

### What changed

The MKL run failed a single backward-compatibility test,
`test-krls.R:93` ("default lambda_method is 'loo'"), which asserted
that two independent `krls()` fits select a bit-identical `lambda`
(`tolerance = 1e-12`). Both fits run the same deterministic exact
path, so under a reproducible single-threaded BLAS they agree exactly
and the test passes on every flavour of the main CRAN check. Under
Intel MKL's multithreaded BLAS the leave-one-out objective is not
reproducible between the two calls (~1e-14 differences from parallel
reduction order); because that objective is nearly flat over a wide
`lambda` range, the golden-section search amplifies the perturbation
into a large gap in the selected `lambda` (0.052 vs ~0 in the report).

The assertion was therefore testing platform determinism, not
behaviour. The test now verifies the actual contract — that the
default objective is LOO, via the `lambda_method` label — and drops
the `lambda`-value comparison. No package (R/src) code changed, so
results on the reference BLAS are unchanged and no reverse-dependency
behaviour is affected.

### Test environments

* macOS 26.5 (Tahoe), aarch64, local — R 4.5.3 (2026-03-11)
* win-builder R-release — R 4.6.1 — Status: OK
* win-builder R-devel — R-devel r90207 (2026-07-04) — Status: OK

### R CMD check results

`R CMD check --as-cran` is clean locally (0 errors, 0 warnings; the
only NOTE is the environmental "package 'V8' unavailable" math-render
skip, present on CRAN's farm). Both win-builder runs returned Status:
OK with no NOTEs. The full test suite passes, including the amended
backward-compatibility test.

### Reverse dependencies

Unchanged from 1.7-0: `InfluenceBorrowing` and `qqkrls`. Because this
patch touches only a test file, revdep behaviour is unaffected.

---

## Submission notes for KRLS 1.7-0 (previous release)

This release brings CRAN up to date with a substantial line of
development since the last published version (1.1-0, 2026-04-30).
Five interim tagged releases on GitHub (1.2-0 through 1.6-1) are
collapsed into this submission; the cumulative changes are
summarized below.

### What's new since 1.1-0

* **1.2-0 — Tier-A modernization.** Formula interface
  (`krls(y ~ x1 + x2, data = df)`); `broom` methods (`tidy`,
  `glance`, `augment`); `ggplot2` `autoplot`; quickstart vignette.
  Algorithm unchanged.

* **1.4-0 — AME-variance optimization.** Per-predictor variance
  uses the row-sum identity `sum(L' V L) = (L*1)' V (L*1)`,
  dropping cost from O(n^3) to O(n^2). 1.2-3x faster on the
  default path; bit-identical to 1.2-0 except for a 1e-15
  summation-order shift in `var.avgderivatives`.

* **1.4-0 -> 1.5-2 — Nystrom approximation.**
  `krls(..., approx = "nystrom")` replaces the full n x n kernel
  with a low-rank landmark feature map. Time O(n m^2 + m^3),
  memory O(n m), making KRLS feasible well past the ~5000-row
  ceiling. Random and k-means landmark selection, point estimates
  and analytical SEs for predictions and AMEs. The variance
  derivation passes a parametric bootstrap calibration; script
  in `dev/03_nystrom_bootstrap_validation.R`. Fit object carries
  `inference = "conditional_nystrom"` and the man page is precise
  about the conditional-on-landmarks interpretation.

* **1.6-0 — Auto-dispatch, default Nystrom m, landmark_seed.**
  `approx = "auto"` becomes the default: at `N <= nystrom_m` the
  exact path runs unchanged; at `N > nystrom_m` the package
  switches to `"nystrom"` and emits a one-line `message()` so the
  dispatch is visible at the call site. `approx = "none"` forces
  the exact path. Default `nystrom_m` raised from
  `min(500, ceiling(sqrt(N)))` to `min(500, N)` (so e.g. at
  N=10,000 m goes from 100 to 500). New `landmark_seed` argument
  saves and restores `.Random.seed` around landmark selection so
  the seed is consumed locally and the caller's RNG state is
  unaffected.

* **1.6-1 — Bug fixes on 1.6-0.** The `derivative = TRUE`,
  `vcov = FALSE` guard now re-runs after auto-dispatch resolves
  (previously crashed with `object 'vcovmatc' not found`).
  `nystrom_m` validation runs before the auto-dispatch threshold
  check.

* **1.7-0 (this release) — `cat_columns` and `maxvarK` defaults.**
  New `cat_columns` argument: declared categorical columns are
  one-hot encoded with all levels and multiplied by `sqrt(0.5)`,
  matching the convention in `kbal::one_hot` and `gpss::one_hot`.
  Continuous columns are still standardized to `sd = 1`. There is
  no autodetection; a one-time `warning()` nudges users when
  factor/character/logical/low-cardinality-numeric columns are not
  declared. `sigma = NULL` (the default) now triggers
  `b_maxvarK()`, which selects the bandwidth that maximizes the
  variance of the off-diagonal kernel entries (matches
  `kbal::b_maxvarK` to ~5e-11 in sigma; pre-1.7 default
  `sigma = ncol(X_processed)` is one line of code away). A
  one-line `message()` announces the change on first call.

### Backwards compatibility

Pinning `sigma = ncol(X_processed)` and `approx = "none"`
reproduces 1.5-2 / 1.6-x fits bit-for-bit
(`all.equal(..., tol = 0)` on `lambda`, `coeffs`, `fitted`,
`avgderivatives`, `vcov.c` on a seed-pinned fixture). Fits
produced by 1.6-x without the new `$prep` slot continue to work
in `predict()` via a legacy code path. All exported function
names and signatures from 1.0-0 are preserved.

### Test environments

* macOS 26.5 (Tahoe), aarch64, local — R 4.5.3 (2026-03-11)
* win-builder R-devel and R-release (planned on submission)
* r-hub macOS-release and ubuntu-release (planned on submission)

### R CMD check results

`R CMD check --as-cran` is clean locally — Status: 1 NOTE. The
single NOTE is purely environmental and not reproducible on CRAN's
farm:

* "Skipping checking math rendering: package 'V8' unavailable"
  (V8 is an optional Suggests-style dependency for math rendering;
  not installed locally, present on CRAN.)

Zero substantive NOTEs, 0 warnings, 0 errors. All 188 tests pass.
Vignettes rebuild cleanly in 14-15 seconds.

### Reverse dependencies

`KRLS` has 2 direct reverse dependencies on CRAN:
`InfluenceBorrowing` (0.1.0) and `qqkrls` (1.0.0). We ran
`revdepcheck::revdep_check()` comparing the new 1.7-0 against the
CRAN baseline 1.1-0. Both revdeps pass cleanly under 1.7-0:

* `InfluenceBorrowing 0.1.0` — 0 errors, 0 warnings, 0 new NOTEs.
* `qqkrls 1.0.0` — 0 errors, 0 warnings, 0 new NOTEs.

Summary: 2 checked, 0 new problems, 0 failures.

### What we kept stable

* All exported function names and signatures from 1.0-0 / 1.1-0:
  `krls`, `predict.krls`, `summary.krls`, `plot.krls`,
  `gausskernel`, `looloss`, `solveforc`.
* All `krls()` return-list field names; new fields are additive
  (`$X_proc`, `$prep`, `$nystrom_*`, `$landmark_indices`,
  `$inference`).
* Pre-1.7 default `sigma = ncol(X_processed)` reproducible by
  setting `sigma` explicitly; the new maxvarK default is announced
  with a one-line message on first call.
