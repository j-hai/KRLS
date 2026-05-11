# =============================================================================
# Nystrom AME standard-error calibration check
# =============================================================================
#
# Under approx = "nystrom" with vcov = TRUE, KRLS reports *conditional
# approximate* standard errors for the average marginal effects: they
# condition on the selected landmarks, fixed lambda, and the low-rank
# feature approximation. This script verifies the analytical SE formula
# is correctly calibrated under that inferential contract via a
# parametric bootstrap that holds those quantities fixed and draws fresh
# Gaussian noise.
#
# Run as a non-CRAN developer validation:
#
#   Rscript dev/03_nystrom_bootstrap_validation.R
#
# Expected behaviour: per-predictor bootstrap SDs of AMEs should match
# the reported var.avgderivatives^0.5 to within Monte Carlo noise
# (roughly 5% for B = 200). A consistent gap of more than ~10% across
# predictors would point to a bug in the variance derivation.

# Validate the source tree rather than any stale installed package.
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(KRLS))
  message("pkgload not available; falling back to installed KRLS ",
          packageVersion("KRLS"))
}

run_check <- function(n = 300, m = 40, d = 3, sigma_noise = 0.3,
                      B = 200, seed = 20260511L) {
  set.seed(seed)
  X <- matrix(rnorm(n * d), n, d)
  colnames(X) <- paste0("x", seq_len(d))
  true_f <- function(X) sin(X[, 1]) + 0.3 * X[, 2] + 0.05 * X[, 3]^2
  y <- true_f(X) + rnorm(n, sd = sigma_noise)

  # Fit once: this locks in landmarks, lambda, and the sigma^2 estimate
  # that the analytical SE conditions on.
  fit <- krls(X, y, approx = "nystrom", nystrom_m = m, print.level = 0)

  reported_se <- sqrt(as.vector(fit$var.avgderivatives))
  fhat        <- as.vector(fit$fitted)
  sigmasq_hat <- mean((y - fhat)^2)

  amebs <- matrix(NA_real_, B, ncol(X))
  for (b in seq_len(B)) {
    y_b   <- fhat + rnorm(n, sd = sqrt(sigmasq_hat))
    fit_b <- krls(X, y_b, approx = "nystrom",
                  landmarks = fit$landmark_indices,
                  lambda    = fit$lambda,
                  print.level = 0)
    amebs[b, ] <- as.vector(fit_b$avgderivatives)
  }
  boot_sd <- apply(amebs, 2, sd)

  data.frame(
    predictor    = colnames(X),
    reported_SE  = round(reported_se, 5),
    bootstrap_SD = round(boot_sd, 5),
    ratio        = round(boot_sd / reported_se, 3)
  )
}

result <- run_check()
cat("\nConditional Nystrom AME SE calibration check\n")
cat("n = 300, m = 40, sigma = 0.3, B = 200 parametric bootstrap draws\n\n")
print(result, row.names = FALSE)
cat("\nReady ratios should sit in roughly [0.9, 1.1] under correct calibration.\n")
