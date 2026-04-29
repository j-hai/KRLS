# =============================================================================
# Capture golden-baseline outputs from KRLS 1.0-0 (frozen reference)
# =============================================================================
#
# Install the frozen 1.0-0 tarball into a temporary library, run a fixed-seed
# set of toy examples covering krls() + summary() + predict() + helpers, and
# save outputs to RDS for later regression testing.
#
# How to run:
#   Rscript dev/01_capture_baseline.R
# =============================================================================

# ---- locate paths -----------------------------------------------------------
script_dir <- if (sys.nframe() > 0) {
  tryCatch(dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
           error = function(e) getwd())
} else getwd()
if (!file.exists(file.path(script_dir, "KRLS_1.0-0_baseline.tar.gz"))) {
  if (file.exists("dev/KRLS_1.0-0_baseline.tar.gz")) {
    script_dir <- file.path(getwd(), "dev")
  }
}
tarball <- file.path(script_dir, "KRLS_1.0-0_baseline.tar.gz")
stopifnot(file.exists(tarball))

baseline_rds <- file.path(script_dir, "baseline_1.0-0.rds")
baseline_log <- file.path(script_dir, "baseline_1.0-0.log")

# ---- install into isolated library -----------------------------------------
templib <- tempfile("krls_baseline_lib_")
dir.create(templib)
cat("Installing", basename(tarball), "into", templib, "\n")
install.packages(tarball, repos = NULL, type = "source", lib = templib,
                 quiet = TRUE, INSTALL_opts = "--no-multiarch")

.libPaths(c(templib, .libPaths()))
suppressPackageStartupMessages(library(KRLS, lib.loc = templib))
stopifnot(packageVersion("KRLS") == "1.0-0")

# ---- deterministic test scenarios ------------------------------------------
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
scenarios <- list()

# Scenario 1: linear additive truth, default settings -----------------------
{
  set.seed(20260429L)
  N  <- 100
  x1 <- rnorm(N)
  x2 <- rbinom(N, 1, 0.3)
  y  <- x1 + 0.5 * x2 + rnorm(N, 0, 0.2)
  X  <- cbind(x1, x2)
  fit <- krls(X = X, y = y, print.level = 0)
  scenarios$s1_linear_default <- list(
    inputs = list(X = X, y = y),
    coeffs        = fit$coeffs,
    fitted        = fit$fitted,
    sigma         = fit$sigma,
    lambda        = fit$lambda,
    R2            = fit$R2,
    avgderivative = if (!is.null(fit$avgderivatives)) fit$avgderivatives else NA,
    derivatives_dim = if (!is.null(fit$derivatives)) dim(fit$derivatives) else NA
  )
}

# Scenario 2: nonlinear truth (x1^3) -----------------------------------------
{
  set.seed(20260429L)
  N  <- 100
  x1 <- rnorm(N)
  x2 <- rbinom(N, 1, 0.3)
  y  <- x1^3 + 0.5 * x2 + rnorm(N, 0, 0.2)
  X  <- cbind(x1, x2)
  fit <- krls(X = X, y = y, print.level = 0)
  scenarios$s2_nonlinear <- list(
    inputs = list(X = X, y = y),
    coeffs        = fit$coeffs,
    fitted        = fit$fitted,
    sigma         = fit$sigma,
    lambda        = fit$lambda,
    R2            = fit$R2,
    avgderivative = if (!is.null(fit$avgderivatives)) fit$avgderivatives else NA
  )
}

# Scenario 3: predict() on new data ------------------------------------------
{
  set.seed(20260429L)
  N  <- 80
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  y  <- x1 + 0.5 * x2 + rnorm(N, 0, 0.2)
  X  <- cbind(x1, x2)
  fit <- krls(X = X, y = y, derivative = FALSE, vcov = FALSE,
              print.level = 0)

  set.seed(987L)
  Xnew <- cbind(rnorm(20), rnorm(20))
  pred <- predict(fit, newdata = Xnew, se.fit = FALSE)
  scenarios$s3_predict <- list(
    fitted    = fit$fitted,
    pred      = pred$fit,
    Xnew      = Xnew
  )
}

# Scenario 4: helpers --------------------------------------------------------
{
  set.seed(20260429L)
  X <- matrix(rnorm(60), 20, 3)
  K <- gausskernel(X = X, sigma = 3)
  scenarios$s4_helpers <- list(
    K_dim    = dim(K),
    K_diag   = diag(K),
    K_first5 = K[1:5, 1:5]
  )
}

# ---- save -------------------------------------------------------------------
saveRDS(scenarios, baseline_rds, version = 2)
sink(baseline_log)
cat("KRLS version:  ", as.character(packageVersion("KRLS")), "\n", sep = "")
cat("R version:     ", R.version.string, "\n", sep = "")
cat("captured at:   ", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
cat("RDS:           ", baseline_rds, "\n", sep = "")
cat("scenarios:     ", length(scenarios), " (", paste(names(scenarios), collapse = ", "), ")\n", sep = "")
cat("\n--- sessionInfo ---\n")
print(sessionInfo())
sink()
cat("Wrote", baseline_rds, "\n")
cat("Wrote", baseline_log, "\n")
