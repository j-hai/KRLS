## Internal helpers for krls(..., approx = "nystrom").
##
## Nystrom replaces the full n x n kernel by a low-rank approximation
## anchored at m << n landmarks. The ridge problem then lives in
## m-dimensional space:
##
##   Phi      = C %*% U %*% diag(D_reg^{-1/2})       (n x m feature map)
##   beta_hat = (Phi' Phi + lambda I)^{-1} Phi' y
##   f_hat(x) = K(x, Z) %*% alpha,   alpha = U %*% (D_reg^{-1/2} * beta)
##
## where C = K(X, Z), W = K(Z, Z), and (U, D) is the eigendecomposition
## of W with relative-ridge stabilization D_reg = pmax(D, eps * max(D)).
##
## See plans/krls-nystrom-1.4-0.md for the full design.

.select_landmarks_random <- function(n, m) {
  sort.int(sample.int(n, m))
}

.validate_nystrom_m <- function(nystrom_m, n) {
  if (!is.numeric(nystrom_m) || length(nystrom_m) != 1L ||
      is.na(nystrom_m) || !is.finite(nystrom_m) ||
      nystrom_m != as.integer(nystrom_m)) {
    stop("nystrom_m must be a finite integer")
  }
  nystrom_m <- as.integer(nystrom_m)
  if (nystrom_m < 1L || nystrom_m > n) {
    stop("nystrom_m must satisfy 1 <= nystrom_m <= nrow(X)")
  }
  nystrom_m
}

.validate_nystrom_eps <- function(nystrom_eps) {
  if (!is.numeric(nystrom_eps) || length(nystrom_eps) != 1L ||
      is.na(nystrom_eps) || !is.finite(nystrom_eps) ||
      nystrom_eps <= 0) {
    stop("nystrom_eps must be a finite positive scalar")
  }
  nystrom_eps
}

## X_std is the already-standardized training matrix; X_centers / X_scales
## are the column means / SDs of the original-scale training X. Index- and
## NULL-form landmarks resolve via X_std rows (already standardized).
## Matrix-form landmarks are taken as user-supplied original-scale
## coordinates and standardized with the same centers/scales so they live
## in the same kernel space the model uses internally.
.resolve_landmarks <- function(landmarks, landmark_method, nystrom_m,
                               X_std, X_centers, X_scales) {
  n <- nrow(X_std)
  d <- ncol(X_std)

  if (is.null(landmarks)) {
    if (is.null(nystrom_m)) {
      nystrom_m <- min(500L, as.integer(ceiling(sqrt(n))))
    }
    nystrom_m <- .validate_nystrom_m(nystrom_m, n)
    if (landmark_method == "kmeans") {
      # k-means on standardized X. Centers come back in the same
      # standardized space, so they slot in alongside index-form
      # landmarks without further transformation. Not bit-stable
      # across R versions even under set.seed() because of the
      # Hartigan-Wong algorithm's initialization; documented.
      km <- stats::kmeans(X_std, centers = nystrom_m, nstart = 10L,
                          iter.max = 50L)
      Z  <- unname(km$centers)
      attr(Z, "kmeans_iter")  <- km$iter
      attr(Z, "kmeans_tot")   <- km$tot.withinss
      return(list(indices = NULL, matrix = Z, method_used = "kmeans"))
    }
    idx <- .select_landmarks_random(n, nystrom_m)
    return(list(indices = idx,
                matrix  = X_std[idx, , drop = FALSE],
                method_used = "random"))
  }

  if (is.numeric(landmarks) && is.null(dim(landmarks))) {
    if (anyNA(landmarks) || any(!is.finite(landmarks)) ||
        any(landmarks != as.integer(landmarks))) {
      stop("landmarks (as indices) must be integer-valued")
    }
    idx <- as.integer(landmarks)
    if (anyNA(idx) || any(idx < 1L) || any(idx > n) || anyDuplicated(idx)) {
      stop("landmarks (as indices) must be unique integers in 1:nrow(X)")
    }
    return(list(indices = idx,
                matrix  = X_std[idx, , drop = FALSE],
                method_used = "user_indices"))
  }

  if (is.matrix(landmarks) || is.data.frame(landmarks)) {
    Z <- as.matrix(landmarks)
    if (!is.numeric(Z))       stop("landmarks matrix must be numeric")
    if (ncol(Z) != d)         stop("ncol(landmarks) must equal ncol(X)")
    if (anyNA(Z) || any(!is.finite(Z)))
      stop("landmarks matrix must contain only finite, non-NA values")
    if (anyDuplicated(Z))
      stop("landmarks matrix must not contain duplicate rows ",
           "(W would be rank-deficient)")
    # User-supplied landmarks are in original X-scale; standardize using
    # the training X's centers/scales to land in the same standardized
    # kernel space the model uses internally.
    Z_std <- scale(Z, center = X_centers, scale = X_scales)
    attr(Z_std, "scaled:center") <- NULL
    attr(Z_std, "scaled:scale")  <- NULL
    return(list(indices = NULL, matrix = Z_std, method_used = "user_matrix"))
  }

  stop("`landmarks` must be NULL, an integer vector of row indices, ",
       "or an m x d numeric matrix (in original X-scale)")
}

## Cross-kernel K(A, B) for a Gaussian kernel with bandwidth sigma.
## Uses the ||a-b||^2 = ||a||^2 + ||b||^2 - 2 a'b decomposition.
.gauss_cross_kernel <- function(A, B, sigma) {
  aa <- rowSums(A * A)
  bb <- rowSums(B * B)
  D2 <- outer(aa, bb, `+`) - 2 * tcrossprod(A, B)
  D2[D2 < 0] <- 0      # guard against tiny negatives from rounding
  exp(-D2 / sigma)
}

## Bound-finding for the Nystrom lambda search. Mirrors lambdasearch()'s
## heuristic but operates on the m-space ridge spectrum (Sigma^2 of Phi),
## not the full-K spectrum.
.nystrom_lambda_bounds <- function(Sigma2, n, L, U) {
  if (is.null(U)) {
    U <- n
    while (sum(Sigma2 / (Sigma2 + U)) < 1) U <- U - 1
  } else {
    stopifnot(is.numeric(U), length(U) == 1L, U > 0)
  }
  if (is.null(L)) {
    q <- which.min(abs(Sigma2 - (max(Sigma2) / 1000)))
    L <- .Machine$double.eps
    while (sum(Sigma2 / (Sigma2 + L)) > q) L <- L + 0.05
  } else {
    stopifnot(is.numeric(L), length(L) == 1L, L >= 0)
  }
  list(L = L, U = U)
}

## Core fit. Operates on already-standardized X and y, returns
## standardized-scale outputs plus the Nystrom bookkeeping that krls()
## will splice into the final fit object.
.fit_krls_nystrom <- function(X, y, sigma, landmarks_resolved,
                              lambda, nystrom_eps,
                              L, U, tol, noisy,
                              compute_vcov = TRUE) {
  n <- nrow(X)
  Z <- landmarks_resolved$matrix
  m <- nrow(Z)
  y_vec <- as.vector(y)
  nystrom_eps <- .validate_nystrom_eps(nystrom_eps)

  C <- .gauss_cross_kernel(X, Z, sigma)            # n x m
  W <- .gauss_cross_kernel(Z, Z, sigma)            # m x m

  We       <- eigen(W, symmetric = TRUE)
  Dvals    <- We$values
  Dvals    <- pmax(Dvals, 0)                       # numerical guard
  Dmax     <- max(Dvals)
  if (Dmax <= 0) stop("anchor kernel W is numerically zero; try a larger sigma")
  D_reg    <- pmax(Dvals, nystrom_eps * Dmax)
  Dinvsqrt <- 1 / sqrt(D_reg)

  # Diagnostics: how many landmark-kernel eigenvalues hit the relative
  # ridge floor, and the pre-floor spectrum range. Useful for spotting
  # near-collinear landmarks / oversized m.
  floored_count <- sum(Dvals < nystrom_eps * Dmax)
  D_min_raw     <- min(Dvals)

  # Phi = C %*% U %*% diag(Dinvsqrt)   (n x m)
  Phi      <- C %*% sweep(We$vectors, 2, Dinvsqrt, `*`)

  # Cache SVD for the lambda loop.
  svd_phi  <- svd(Phi)
  Sigma2   <- svd_phi$d^2

  ## LOO loss at a given lambda, O(n m) using the cached SVD.
  loo_loss <- function(lambda_val) {
    w     <- Sigma2 / (Sigma2 + lambda_val)
    yfit  <- as.numeric(svd_phi$u %*% (w * crossprod(svd_phi$u, y_vec)))
    diagS <- as.numeric((svd_phi$u^2) %*% w)
    # Guard against (1 - diagS) ~ 0 from rounding when lambda is tiny
    # and the leave-one-out residual blows up; replace with Inf so the
    # search avoids that region.
    denom <- 1 - diagS
    if (any(denom <= .Machine$double.eps)) return(.Machine$double.xmax)
    sum(((y_vec - yfit) / denom)^2)
  }

  if (is.null(lambda)) {
    if (is.null(tol)) tol <- 1e-3 * n
    bounds <- .nystrom_lambda_bounds(Sigma2, n, L, U)
    L <- bounds$L; U <- bounds$U

    # Golden-section search, structurally identical to lambdasearch().
    gr  <- 0.381966
    X1  <- L + gr * (U - L); X2 <- U - gr * (U - L)
    S1  <- loo_loss(X1);     S2 <- loo_loss(X2)
    if (noisy) cat("L:", L, "X1:", X1, "X2:", X2, "U:", U,
                   "S1:", S1, "S2:", S2, "\n")
    while (abs(S1 - S2) > tol) {
      if (S1 < S2) {
        U  <- X2; X2 <- X1; X1 <- L + gr * (U - L)
        S2 <- S1; S1 <- loo_loss(X1)
      } else {
        L  <- X1; X1 <- X2; X2 <- U - gr * (U - L)
        S1 <- S2; S2 <- loo_loss(X2)
      }
      if (noisy) cat("L:", L, "X1:", X1, "X2:", X2, "U:", U,
                     "S1:", S1, "S2:", S2, "\n")
    }
    lambda <- if (S1 < S2) X1 else X2
    if (noisy) cat("Lambda:", lambda, "\n")
  } else {
    stopifnot(is.numeric(lambda), length(lambda) == 1L, lambda > 0)
  }

  # Solution at chosen lambda.
  w        <- Sigma2 / (Sigma2 + lambda)
  yfit_std <- as.numeric(svd_phi$u %*% (w * crossprod(svd_phi$u, y_vec)))
  # beta = V (Sigma / (Sigma^2 + lambda)) U' y     (length m)
  beta     <- as.numeric(
    svd_phi$v %*% ((svd_phi$d * crossprod(svd_phi$u, y_vec)) / (Sigma2 + lambda))
  )
  # alpha = U diag(Dinvsqrt) beta    (length m; gives f(x) = K(x,Z) %*% alpha)
  alpha    <- as.numeric(We$vectors %*% (Dinvsqrt * beta))

  diagS    <- as.numeric((svd_phi$u^2) %*% w)
  Looe_std <- sum(((y_vec - yfit_std) / (1 - diagS))^2)
  vcov_alpha_std <- NULL
  if (isTRUE(compute_vcov)) {
    sigmasq_std <- as.numeric((1 / n) * crossprod(y_vec - yfit_std))
    beta_weights <- sigmasq_std * Sigma2 / (Sigma2 + lambda)^2
    vcov_beta <- tcrossprod(sweep(svd_phi$v, 2, beta_weights, `*`),
                            svd_phi$v)
    alpha_map <- sweep(We$vectors, 2, Dinvsqrt, `*`)
    vcov_alpha_std <- tcrossprod(alpha_map %*% vcov_beta, alpha_map)
  }

  list(
    coeffs           = matrix(alpha, ncol = 1),   # length-m
    fitted_std       = yfit_std,                  # length-n
    lambda           = lambda,
    Looe_std         = Looe_std,
    landmarks        = Z,                         # m x d, standardized
    landmark_indices = landmarks_resolved$indices,
    landmark_method  = landmarks_resolved$method_used,
    W_eigen          = list(values = Dvals, vectors = We$vectors),
    Dinvsqrt         = Dinvsqrt,
    Sigma2           = Sigma2,                    # m-length, for eff_df / GCV
    vcov_alpha_std   = vcov_alpha_std,
    C                = C,                         # n x m, for derivative loop
    nystrom_m        = m,
    nystrom_eps      = nystrom_eps,
    floored_count    = floored_count,
    D_min_raw        = D_min_raw,
    D_max_raw        = Dmax
  )
}

## Pointwise derivatives under Nystrom:
##   d/dx_d f_hat(x_i) = sum_j alpha_j (-2/sigma) (X[i,d] - Z[j,d]) C[i,j]
.nystrom_derivatives <- function(X, Z, C, alpha, sigma) {
  n <- nrow(X); d <- ncol(X); m <- nrow(Z)
  derivmat <- matrix(NA_real_, n, d)
  scale_factor <- -2 / sigma
  for (k in seq_len(d)) {
    # distk[i, j] = X[i, k] - Z[j, k]
    distk <- outer(X[, k], Z[, k], `-`)         # n x m
    Lmat  <- distk * C                          # n x m
    derivmat[, k] <- scale_factor * (Lmat %*% alpha)
  }
  derivmat
}

.nystrom_ame_variances <- function(X, Z, C, vcov_alpha, sigma) {
  if (is.null(vcov_alpha)) return(NULL)
  d <- ncol(X)
  out <- matrix(NA_real_, 1, d)
  scale_factor <- -2 / sigma
  for (k in seq_len(d)) {
    distk <- outer(X[, k], Z[, k], `-`)
    Lmat <- distk * C
    g <- as.numeric(scale_factor * colMeans(Lmat))
    out[1, k] <- as.numeric(crossprod(g, vcov_alpha %*% g))
  }
  out
}
