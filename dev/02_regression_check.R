# =============================================================================
# Regression check: current dev source vs. frozen KRLS 1.0-0 baseline
# =============================================================================
# Same contract as the ebal/Synth check scripts. Default tol 1e-8.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
opt_tol <- 1e-8
opt_verbose <- FALSE
for (a in args) {
  if (grepl("^--tol=", a)) opt_tol <- as.numeric(sub("^--tol=", "", a))
  if (a == "--verbose") opt_verbose <- TRUE
}

expected_diffs <- list(
  # Filled in as Phase 2 fixes land.
)

pkg_root <- getwd()
stopifnot(file.exists(file.path(pkg_root, "DESCRIPTION")))
baseline_rds <- file.path(pkg_root, "dev", "baseline_1.0-0.rds")
stopifnot(file.exists(baseline_rds))

templib <- tempfile("krls_dev_lib_")
dir.create(templib)
cat("Installing current source from", pkg_root, "into temp lib...\n")
res <- system2(file.path(R.home("bin"), "R"),
               c("CMD", "INSTALL", "--no-multiarch", "-l", shQuote(templib),
                 shQuote(pkg_root)),
               stdout = TRUE, stderr = TRUE)
if (!is.null(attr(res, "status")) && attr(res, "status") != 0) {
  cat(res, sep = "\n"); stop("R CMD INSTALL failed")
}

.libPaths(c(templib, .libPaths()))
suppressPackageStartupMessages(library(KRLS, lib.loc = templib))
cat("Loaded KRLS", as.character(packageVersion("KRLS")), "from temp lib\n\n")

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

run_scenarios <- function() {
  out <- list()

  # s1
  set.seed(20260429L)
  N  <- 100
  x1 <- rnorm(N); x2 <- rbinom(N, 1, 0.3)
  y  <- x1 + 0.5 * x2 + rnorm(N, 0, 0.2)
  X  <- cbind(x1, x2)
  fit <- krls(X = X, y = y, print.level = 0)
  out$s1_linear_default <- list(
    inputs = list(X = X, y = y),
    coeffs        = fit$coeffs,
    fitted        = fit$fitted,
    sigma         = fit$sigma,
    lambda        = fit$lambda,
    R2            = fit$R2,
    avgderivative = if (!is.null(fit$avgderivatives)) fit$avgderivatives else NA,
    derivatives_dim = if (!is.null(fit$derivatives)) dim(fit$derivatives) else NA
  )

  # s2
  set.seed(20260429L)
  N  <- 100
  x1 <- rnorm(N); x2 <- rbinom(N, 1, 0.3)
  y  <- x1^3 + 0.5 * x2 + rnorm(N, 0, 0.2)
  X  <- cbind(x1, x2)
  fit <- krls(X = X, y = y, print.level = 0)
  out$s2_nonlinear <- list(
    inputs = list(X = X, y = y),
    coeffs        = fit$coeffs,
    fitted        = fit$fitted,
    sigma         = fit$sigma,
    lambda        = fit$lambda,
    R2            = fit$R2,
    avgderivative = if (!is.null(fit$avgderivatives)) fit$avgderivatives else NA
  )

  # s3
  set.seed(20260429L)
  N  <- 80
  x1 <- rnorm(N); x2 <- rnorm(N)
  y  <- x1 + 0.5 * x2 + rnorm(N, 0, 0.2)
  X  <- cbind(x1, x2)
  fit <- krls(X = X, y = y, derivative = FALSE, vcov = FALSE,
              print.level = 0)
  set.seed(987L)
  Xnew <- cbind(rnorm(20), rnorm(20))
  pred <- predict(fit, newdata = Xnew, se.fit = FALSE)
  out$s3_predict <- list(
    fitted    = fit$fitted,
    pred      = pred$fit,
    Xnew      = Xnew
  )

  # s4
  set.seed(20260429L)
  X <- matrix(rnorm(60), 20, 3)
  K <- gausskernel(X = X, sigma = 3)
  out$s4_helpers <- list(
    K_dim    = dim(K),
    K_diag   = diag(K),
    K_first5 = K[1:5, 1:5]
  )

  out
}

current  <- suppressWarnings(suppressMessages(run_scenarios()))
baseline <- readRDS(baseline_rds)

compare_one <- function(cur, base, path = "", tol = opt_tol) {
  if (is.list(cur) && is.list(base)) {
    diffs <- character()
    for (k in names(base)) {
      diffs <- c(diffs, compare_one(cur[[k]], base[[k]],
                                    paste0(path, "$", k), tol))
    }
    return(diffs)
  }
  if (is.numeric(cur) && is.numeric(base)) {
    if (length(cur) != length(base)) {
      return(sprintf("%s: length differs (cur=%d base=%d)", path,
                     length(cur), length(base)))
    }
    if (length(cur) == 0L) return(character())
    md <- max(abs(cur - base), na.rm = TRUE)
    if (anyNA(cur) != anyNA(base)) {
      return(sprintf("%s: NA pattern differs", path))
    }
    if (is.finite(md) && md > tol) {
      return(sprintf("%s: max |diff| = %.3g > tol=%.1g", path, md, tol))
    }
    return(character())
  }
  if (!identical(cur, base)) return(sprintf("%s: identical() failed", path))
  character()
}

cat("Tolerance:", opt_tol, "\n")
cat("Comparing", length(current), "scenarios:\n\n")

unexpected <- 0L; expected_changed <- 0L; ok <- 0L
for (nm in names(baseline)) {
  if (!nm %in% names(current)) {
    cat(sprintf("  [MISSING] %s\n", nm)); unexpected <- unexpected + 1L; next
  }
  diffs <- compare_one(current[[nm]], baseline[[nm]], path = nm)
  if (length(diffs) == 0L) {
    cat(sprintf("  [OK]      %s\n", nm)); ok <- ok + 1L
  } else if (nm %in% names(expected_diffs)) {
    cat(sprintf("  [EXPECTED] %s — %s\n", nm, expected_diffs[[nm]]))
    if (opt_verbose) for (d in diffs) cat(sprintf("            %s\n", d))
    expected_changed <- expected_changed + 1L
  } else {
    cat(sprintf("  [DIFF]    %s\n", nm))
    for (d in diffs) cat(sprintf("            %s\n", d))
    unexpected <- unexpected + 1L
  }
}

cat(sprintf("\nSummary: %d ok, %d expected-changed, %d unexpected\n",
            ok, expected_changed, unexpected))
if (unexpected > 0L) quit(save = "no", status = 1L)
quit(save = "no", status = 0L)
