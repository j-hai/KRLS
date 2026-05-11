## Tests for the Nyström approximation path (approx = "nystrom").
##
## Scope: predictions, fitted values, AMEs, conditional approximate
## inference, and the basic API contract.

test_that("Nyström fit object has the expected structure", {
  set.seed(1)
  d <- make_nonlinear_data(N = 100)
  fit <- krls(d$X, d$y, approx = "nystrom", nystrom_m = 30, print.level = 0)

  expect_s3_class(fit, "krls")
  expect_equal(fit$approx, "nystrom")
  expect_equal(fit$nystrom_m, 30L)
  expect_equal(nrow(fit$landmarks), 30L)
  expect_equal(ncol(fit$landmarks), ncol(d$X))
  expect_length(as.vector(fit$coeffs), 30L)

  # Conditional approximate inference is available under vcov=TRUE,
  # while the full n x n fitted-value covariance is not stored.
  expect_equal(dim(fit$vcov.c), c(30L, 30L))
  expect_null(fit$vcov.fitted)
  expect_equal(dim(fit$var.avgderivatives), c(1L, ncol(d$X)))
  expect_true(all(is.finite(fit$var.avgderivatives)))
  expect_equal(fit$inference, "conditional_nystrom")
  expect_null(fit$K)            # full kernel not stored (memory win)

  # Point estimates present.
  expect_equal(dim(fit$derivatives), c(nrow(d$X), ncol(d$X)))
  expect_equal(dim(fit$avgderivatives), c(1, ncol(d$X)))
})

test_that("predict() works under Nyström, including conditional SEs", {
  set.seed(2)
  d <- make_nonlinear_data(N = 100)
  fit <- krls(d$X, d$y, approx = "nystrom", nystrom_m = 40, print.level = 0)

  newx <- d$X[1:5, , drop = FALSE]
  p <- predict(fit, newdata = newx)
  expect_length(as.vector(p$fit), 5L)
  expect_null(p$se.fit)
  expect_null(p$vcov.fit)

  p_se <- predict(fit, newdata = newx, se.fit = TRUE)
  expect_length(as.vector(p_se$se.fit), 5L)
  expect_equal(dim(p_se$vcov.fit), c(5L, 5L))
  expect_true(all(is.finite(p_se$se.fit)))
  expect_true(all(as.vector(p_se$se.fit) >= 0))

  fit_no_v <- krls(d$X, d$y, approx = "nystrom", nystrom_m = 40,
                   vcov = FALSE, print.level = 0)
  expect_error(
    predict(fit_no_v, newdata = newx, se.fit = TRUE),
    "vcov = TRUE"
  )
})

test_that("Nyström with m = n produces a quality fit comparable to exact", {
  # m = n is NOT bit-equivalent to exact KRLS: the relative-ridge
  # stabilization (D_reg = pmax(D, eps * max(D))) floors numerical-
  # noise eigenvalues that exact KRLS can exploit when picking very
  # small lambda. Both are sensible fits; we test fit quality, not
  # bit-parity. (Users who want near-bit-parity can pass nystrom_eps
  # near .Machine$double.eps.)
  d <- make_nonlinear_data(N = 60)

  fit_e <- krls(d$X, d$y, print.level = 0)
  fit_n <- krls(d$X, d$y, approx = "nystrom",
                landmarks = seq_len(nrow(d$X)),
                print.level = 0)

  expect_gt(fit_n$R2, 0.95)
  # Fits agree in shape if not in exact values.
  expect_gt(cor(as.vector(fit_e$fitted), as.vector(fit_n$fitted)), 0.99)
  # AMEs in the same ballpark.
  expect_gt(cor(as.vector(fit_e$avgderivatives),
                as.vector(fit_n$avgderivatives)), 0.95)
})

test_that("Nyström runs and predicts cleanly at m << n", {
  set.seed(4)
  d <- make_nonlinear_data(N = 300)
  fit <- krls(d$X, d$y, approx = "nystrom", nystrom_m = 30, print.level = 0)

  expect_gt(fit$R2, 0.9)          # should still fit reasonably
  yhat <- predict(fit, newdata = d$X)$fit
  expect_gt(cor(as.vector(yhat), d$y), 0.95)
})

test_that("landmarks accepts NULL, integer indices, and an m x d matrix", {
  set.seed(5)
  d <- make_nonlinear_data(N = 80)

  # NULL: auto-select
  fit_auto <- krls(d$X, d$y, approx = "nystrom",
                   nystrom_m = 20, print.level = 0)
  expect_equal(fit_auto$nystrom_m, 20L)
  expect_false(is.null(fit_auto$landmark_indices))

  # Integer indices
  idx <- c(1L, 5L, 10L, 20L, 40L)
  fit_idx <- krls(d$X, d$y, approx = "nystrom",
                  landmarks = idx, print.level = 0)
  expect_equal(fit_idx$nystrom_m, length(idx))
  expect_equal(fit_idx$landmark_indices, idx)
  # The standardized landmarks should equal the standardized rows.
  Xstd <- scale(d$X, center = TRUE, scale = apply(d$X, 2, sd))
  expect_equal(unname(fit_idx$landmarks), unname(Xstd[idx, , drop = FALSE]),
               tolerance = 1e-12)

  # m x d matrix in original X-scale (the natural user input).
  # Indices and an explicit matrix of the same rows must produce the
  # same fit; that was a 1.4-0 regression where matrix landmarks were
  # treated as already-standardized.
  Z_orig <- d$X[c(2, 7, 15), , drop = FALSE]
  fit_mat <- krls(d$X, d$y, approx = "nystrom",
                  landmarks = Z_orig, print.level = 0)
  expect_equal(fit_mat$nystrom_m, 3L)
  expect_null(fit_mat$landmark_indices)
  # Standardized stored landmarks should match the standardized rows
  # selected via indices.
  expect_equal(unname(fit_mat$landmarks),
               unname(Xstd[c(2, 7, 15), , drop = FALSE]),
               tolerance = 1e-12)

  fit_idx_same <- krls(d$X, d$y, approx = "nystrom",
                       landmarks = c(2L, 7L, 15L), print.level = 0)
  expect_equal(as.vector(fit_mat$fitted),
               as.vector(fit_idx_same$fitted),
               tolerance = 1e-10)
})

test_that("summary() and tidy() work on Nystrom fits with and without SEs", {
  d <- make_nonlinear_data(N = 80)
  fit <- krls(d$X, d$y, approx = "nystrom", nystrom_m = 25, print.level = 0)

  # summary() must report conditional approximate inference by default.
  s <- capture.output(summary_obj <- summary(fit))
  expect_true(any(grepl("Average Marginal Effects", s)))
  expect_true("Est" %in% colnames(summary_obj$coefficients))
  expect_true("Std. Error" %in% colnames(summary_obj$coefficients))
  expect_true(all(is.finite(summary_obj$coefficients[, "Std. Error"])))

  # tidy() must produce finite SE columns when vcov=TRUE.
  td <- generics::tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), ncol(d$X))
  expect_true(all(is.finite(td$std.error)))
  expect_true(all(is.finite(td$p.value)))
  expect_true(all(!is.na(td$estimate)))

  # With vcov=FALSE, summaries remain point-estimate-only and tidy()
  # keeps stable inference columns filled with NA.
  fit_no_v <- krls(d$X, d$y, approx = "nystrom", nystrom_m = 25,
                   vcov = FALSE, print.level = 0)
  s_no_v <- capture.output(summary_no_v <- summary(fit_no_v))
  expect_true(any(grepl("Average Marginal Effects", s_no_v)))
  expect_true("Est" %in% colnames(summary_no_v$coefficients))
  expect_false("Std. Error" %in% colnames(summary_no_v$coefficients))
  td_no_v <- generics::tidy(fit_no_v)
  expect_true(all(is.na(td_no_v$std.error)))
  expect_true(all(is.na(td_no_v$p.value)))
})

test_that("Nyström binary first differences get approximate variances", {
  d <- make_linear_data(N = 80)
  fit <- krls(d$X, d$y, approx = "nystrom", nystrom_m = 25, print.level = 0)
  expect_true(fit$binaryindicator[1, "x2"])
  expect_false(is.null(fit$var.avgderivatives))
  expect_true(is.finite(fit$var.avgderivatives[1, "x2"]))
  expect_gt(fit$var.avgderivatives[1, "x2"], 0)
})

test_that("random landmarks are reproducible under set.seed()", {
  d <- make_nonlinear_data(N = 80)
  set.seed(42); fit1 <- krls(d$X, d$y, approx = "nystrom",
                              nystrom_m = 25, print.level = 0)
  set.seed(42); fit2 <- krls(d$X, d$y, approx = "nystrom",
                              nystrom_m = 25, print.level = 0)
  expect_equal(fit1$landmark_indices, fit2$landmark_indices)
  expect_equal(as.vector(fit1$coeffs), as.vector(fit2$coeffs))
})

test_that("Nyström validates malformed landmarks", {
  d <- make_nonlinear_data(N = 50)

  expect_error(
    krls(d$X, d$y, approx = "nystrom", landmarks = c(1, 1, 2), print.level = 0),
    "unique"
  )
  expect_error(
    krls(d$X, d$y, approx = "nystrom", landmarks = c(0L, 1L), print.level = 0),
    "1:nrow"
  )
  expect_error(
    krls(d$X, d$y, approx = "nystrom",
         landmarks = matrix(0, 3, ncol(d$X) + 1), print.level = 0),
    "ncol"
  )
  expect_error(
    krls(d$X, d$y, approx = "nystrom", nystrom_m = 1.5, print.level = 0),
    "finite integer"
  )
  expect_error(
    krls(d$X, d$y, approx = "nystrom", nystrom_m = NA, print.level = 0),
    "finite integer"
  )
  expect_error(
    krls(d$X, d$y, approx = "nystrom", nystrom_eps = 0, print.level = 0),
    "finite positive scalar"
  )
})

test_that("Nyström rejects unsupported options cleanly", {
  d <- make_nonlinear_data(N = 50)

  expect_error(
    krls(d$X, d$y, approx = "nystrom", whichkernel = "linear", print.level = 0),
    "gaussian"
  )
  expect_error(
    krls(d$X, d$y, approx = "nystrom", eigtrunc = 0.001, print.level = 0),
    "eigtrunc"
  )
})

# --- 1.4-1 polish bundle -----------------------------------------------------

test_that("landmark_method = 'kmeans' produces sensible landmarks", {
  d <- make_nonlinear_data(N = 120)
  set.seed(2026)
  fit <- krls(d$X, d$y, approx = "nystrom", nystrom_m = 25,
              landmark_method = "kmeans", print.level = 0)
  expect_equal(fit$landmark_method, "kmeans")
  expect_null(fit$landmark_indices)
  expect_equal(dim(fit$landmarks), c(25L, ncol(d$X)))
  expect_gt(fit$R2, 0.85)
})

test_that("get_landmarks() returns coordinates in the requested scale", {
  d <- make_nonlinear_data(N = 80)
  fit <- krls(d$X, d$y, approx = "nystrom", nystrom_m = 20, print.level = 0)

  Z_orig <- get_landmarks(fit)
  Z_std  <- get_landmarks(fit, scale = "standardized")
  expect_equal(dim(Z_orig), c(20L, ncol(d$X)))
  expect_equal(unname(Z_std), unname(fit$landmarks), tolerance = 1e-12)
  # Round-trip: passing original-scale landmarks back reproduces the fit.
  fit_back <- krls(d$X, d$y, approx = "nystrom",
                   landmarks = Z_orig, print.level = 0)
  expect_equal(as.vector(fit$fitted), as.vector(fit_back$fitted),
               tolerance = 1e-8)
})

test_that("get_landmarks() errors on non-Nyström fits", {
  d <- make_nonlinear_data(N = 60)
  fit_exact <- krls(d$X, d$y, print.level = 0)
  expect_error(get_landmarks(fit_exact), "approx = 'nystrom'")
  expect_error(get_landmarks(list(foo = 1)), "class 'krls'")
})

test_that("glance() Nyström columns and eff_df are populated", {
  d <- make_nonlinear_data(N = 100)
  fit <- krls(d$X, d$y, approx = "nystrom", nystrom_m = 30, print.level = 0)
  g <- generics::glance(fit)
  expect_true(all(c("approx", "nystrom_m", "inference") %in% colnames(g)))
  expect_equal(g$approx, "nystrom")
  expect_equal(g$nystrom_m, 30L)
  expect_equal(g$inference, "conditional_nystrom")
  expect_true(is.finite(g$eff_df))
  expect_true(g$eff_df > 0 && g$eff_df <= 30)

  # Exact path still reports approx = "none" and a finite eff_df.
  fit_e <- krls(d$X, d$y, print.level = 0)
  g_e <- generics::glance(fit_e)
  expect_equal(g_e$approx, "none")
  expect_true(is.na(g_e$nystrom_m))
  expect_true(is.finite(g_e$eff_df))
})

test_that("summary() prints Nyström approximation block with diagnostics", {
  d <- make_nonlinear_data(N = 100)
  fit <- krls(d$X, d$y, approx = "nystrom", nystrom_m = 25, print.level = 0)
  out <- capture.output(summary(fit))
  expect_true(any(grepl("Approximation: Nystrom with m = 25", out)))
  expect_true(any(grepl("Inference: conditional approximate", out)))
  expect_true(any(grepl("Landmark kernel: condition", out)))
})

test_that("L >= U lambda-window is rejected", {
  d <- make_nonlinear_data(N = 50)
  expect_error(krls(d$X, d$y, L = 10, U = 5, print.level = 0),
               "L must be strictly less than U")
  expect_error(krls(d$X, d$y, L = 1, U = 1, print.level = 0),
               "L must be strictly less than U")
})

test_that("malformed landmark matrices are rejected", {
  d <- make_nonlinear_data(N = 60)
  # Duplicate rows make the landmark kernel rank-deficient.
  Z_dup <- rbind(d$X[1, , drop = FALSE], d$X[1, , drop = FALSE],
                 d$X[2, , drop = FALSE])
  expect_error(krls(d$X, d$y, approx = "nystrom", landmarks = Z_dup,
                    print.level = 0),
               "duplicate rows")
  # NA / non-finite entries.
  Z_bad <- d$X[1:3, , drop = FALSE]; Z_bad[1, 1] <- NA
  expect_error(krls(d$X, d$y, approx = "nystrom", landmarks = Z_bad,
                    print.level = 0),
               "finite")
})
