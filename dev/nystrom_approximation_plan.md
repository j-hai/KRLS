# Nyström Approximation Plan

## Goal

Add an explicitly approximate KRLS fit path for larger data. The exact KRLS
implementation remains the default. The Nyström path trades exact kernel
inference for lower memory use and faster prediction-oriented fits.

## API

Proposed arguments to `krls()`:

```r
approx = c("none", "nystrom")
nystrom_m = NULL
landmarks = NULL
landmark_method = c("random", "kmeans")
nystrom_eps = 1e-8
```

`approx` should be user-facing because this is an approximation, not just a
backend. `method` can remain reserved for exact numerical implementations if it
is kept at all.

`landmarks` accepts:

- `NULL`: choose landmarks using `landmark_method`.
- Integer row indices into `X`: common case for sampled training landmarks.
- An `m x d` numeric matrix: supports external anchors, grids, or transferred
  landmark sets.

Store the resolved landmark matrix in the fitted object. If landmarks came from
training rows, also store the integer indices.

## Feature Map

Let `Z` be the `m x d` landmark matrix, `C = K(X, Z)`, and `W = K(Z, Z)`.
Compute:

```text
W = U D U'
D_reg = pmax(D, nystrom_eps * max(D))
Phi = C U diag(D_reg^(-1/2))
```

The relative ridge on `W` is more scale-robust than adding a fixed epsilon.
For `m = n` and landmarks equal to the training rows, this should closely
recover the exact KRLS fit, modulo numerical stabilization.

Fit standardized `y` with ordinary ridge in feature space:

```text
beta(lambda) = solve(Phi' Phi + lambda I, Phi' y)
fitted = Phi beta(lambda)
```

Do not store the full exact-style `n x n` kernel matrix in Nyström fits. Store
only the pieces needed for prediction and diagnostics.

## Lambda Search

The exact `lambdasearch()` bounds are based on the full kernel spectrum. The
Nyström ridge problem lives in `m`-space, so its relevant spectrum is that of
`crossprod(Phi)`.

Either generalize the bound-finding logic to accept arbitrary eigenvalues, or
factor the bound calculation into a helper that both the exact and Nyström paths
can call with the appropriate spectrum.

Use a one-time SVD of `Phi` to make each lambda evaluation cheap:

```text
Phi = U_p Sigma V_p'
Phi' Phi + lambda I = V_p (Sigma^2 + lambda I) V_p'
```

Then:

```text
fitted(lambda) = U_p diag(Sigma^2 / (Sigma^2 + lambda)) U_p' y
diag(S(lambda)) = rowSums(U_p^2 * rep_each(Sigma^2 / (Sigma^2 + lambda)))
loo_resid = residual / (1 - diag(S))
loo_loss = sum(loo_resid^2)
```

With cached `U_p`, `Sigma`, and `U_p' y`, golden-section lambda search becomes
approximately `O(nm)` per evaluation rather than repeatedly solving an
`m x m` system.

GCV can be computed alongside LOO:

```text
tr(S) = sum(Sigma^2 / (Sigma^2 + lambda))
gcv = mean(residual^2) / (1 - tr(S) / n)^2
```

GCV is useful as a fallback if LOO becomes numerically unstable for pathological
inputs.

## Prediction

For new data `X_new`, compute:

```text
C_new = K(X_new, Z)
Phi_new = C_new U diag(D_reg^(-1/2))
pred = Phi_new beta
```

`predict.krls()` should branch on `object$approx`. For Nyström fits,
`se.fit = TRUE` should error clearly until approximate fitted-value covariance
support is deliberately derived and tested.

## Derivatives And Inference

Pointwise derivatives are feasible for Gaussian kernels because:

```text
f(x) = phi(x)' beta
```

Differentiate the landmark kernel row, then multiply by
`U diag(D_reg^(-1/2)) beta`.

First implementation should keep inference conservative:

- `vcov.c = NULL`
- `vcov.fitted = NULL`
- `var.avgderivatives = NULL`
- `summary.krls()` and `tidy.krls()` omit SE/p-value columns when variances are
  unavailable.

Users who need exact derivative SEs should use `approx = "none"`.

## Implementation Phases

1. Add argument validation and fitted-object metadata for `approx`.
2. Implement random training-row landmarks.
3. Implement landmark indices and landmark matrix inputs.
4. Build Gaussian cross-kernel helpers for `K(X, Z)`.
5. Implement the stabilized `W` eigensystem and `Phi` construction.
6. Implement cached-SVD lambda search for Nyström fits.
7. Fit ridge coefficients and fitted values on standardized data.
8. Add `predict.krls()` support for Nyström objects.
9. Add point-derivative support without derivative variances.
10. Add `kmeans` landmarks after the random-landmark path is stable.

## Tests

Core tests:

- `approx = "none"` remains default and exact behavior is unchanged.
- Random landmark fits run and predict without storing an `n x n` kernel.
- Integer landmark indices and numeric landmark matrices both work.
- With `m = n` and landmarks equal to all training rows, Nyström fitted values
  closely match exact KRLS.
- Cached-SVD LOO lambda search matches a direct ridge/hat-diagonal calculation
  on small problems.
- `predict(..., se.fit = TRUE)` errors clearly for Nyström fits.
- Summary and broom methods gracefully omit inference columns when variance
  fields are `NULL`.

Performance tests should compare memory and runtime for moderate-to-large `n`
with fixed `m << n`, but avoid brittle wall-clock assertions in CRAN tests.

