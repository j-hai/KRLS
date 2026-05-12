krls <-
function(     X=NULL,
              y=NULL,
              whichkernel="gaussian",
              lambda=NULL,
              sigma=NULL,
              derivative=TRUE,
              binary=TRUE,
              vcov=TRUE,
              print.level=1,
              L=NULL,
              U=NULL,
              tol=NULL,
							eigtrunc=NULL,
              data=NULL,
              approx=c("auto","none","nystrom"),
              nystrom_m=NULL,
              landmarks=NULL,
              landmark_method=c("random","kmeans"),
              nystrom_eps=sqrt(.Machine$double.eps),
              landmark_seed=NULL,
              lambda_method=c("loo","gcv")){

      approx          <- match.arg(approx)
      landmark_method <- match.arg(landmark_method)
      lambda_method   <- match.arg(lambda_method)

      # ---- formula interface --------------------------------------------------
      # If the user passed a two-sided formula as the first argument, build
      # X and y from formula + data and continue with the matrix interface.
      # Implemented as in-function dispatch (rather than UseMethod) to avoid
      # CRAN R CMD check NOTEs about long-standing top-level functions
      # whose names accidentally match the <generic>.<class> pattern.
      if (inherits(X, "formula")) {
        formula <- X
        if (is.null(data))
          stop("'data' is required when calling krls() with a formula")
        if (length(formula) != 3L)
          stop("formula must be two-sided, e.g. y ~ x1 + x2")
        mf <- model.frame(formula, data = data, na.action = stats::na.pass)
        yvec <- stats::model.response(mf)
        if (is.null(yvec))
          stop("formula has no response (left-hand side); expected y ~ x1 + x2 + ...")
        if (any(is.na(yvec)))
          stop("y contains missing data")
        if (any(is.na(mf)))
          stop("X contains missing data")
        Xmat <- stats::model.matrix(formula, data = mf)
        if ("(Intercept)" %in% colnames(Xmat))
          Xmat <- Xmat[, colnames(Xmat) != "(Intercept)", drop = FALSE]
        X <- Xmat
        y <- yvec
      }

      # checks
      y <- as.matrix(y)
      X <- as.matrix(X)

      if (is.numeric(X)==FALSE){
       stop("X must be numeric")
      }
      if (is.numeric(y)==FALSE){
       stop("y must be numeric")
      }
      if (sum(is.na(X)) > 0) {
        stop("X contains missing data")
      }
      if (sum(is.na(y)) > 0) {
        stop("y contains missing data")
      }
      if (var(y) == 0) {
        stop("y is a constant (does not vary)")
      }

      if (!is.null(eigtrunc)){
          if (!is.numeric(eigtrunc)) stop("eigtrunc, if used, must be numeric")
          if (eigtrunc>1 | eigtrunc<0) stop("eigtrunc must be between 0 and 1")
          if (eigtrunc==0) {
            eigtrunc=NULL
            warning("eigtrunc of 0 equivalent to no eigen truncation")}
      }

      n <- nrow(X)
      d <- ncol(X)
      if (n!=nrow(y)){
       stop("nrow(X) not equal to number of elements in y")
      }

      stopifnot(
                is.logical(derivative),
                is.logical(vcov),
                is.logical(binary)
                )

      # (The derivative=TRUE / vcov=FALSE guard moved past the
      # auto-dispatch block below: the constraint applies to the exact
      # path but not to Nystrom, which returns point-estimate
      # derivatives without SEs. We need to know which path approx
      # will actually resolve to before enforcing the rule.)

      # If the user supplied both L and U for the lambda search, enforce
      # L < U up front; otherwise lambdasearch()'s golden-section step
      # produces nonsense bounds.
      if (!is.null(L) && !is.null(U)) {
        if (!is.numeric(L) || !is.numeric(U) || length(L) != 1L ||
            length(U) != 1L) {
          stop("L and U must each be a single numeric scalar")
        }
        if (!(L < U)) {
          stop("L must be strictly less than U for the lambda search window")
        }
      }


      # default sigma to dim of X
      if(is.null(sigma)) { sigma <- d
      } else {
        stopifnot(is.vector(sigma),
                  length(sigma)==1,
                  is.numeric(sigma),
                  sigma>0)
      }

      # column names
      if (is.null(colnames(X))) {
      colnames(X) <- paste("x",1:d,sep="")
      }

      # scale
      X.init <- X
      X.init.sd <- apply(X.init,2,sd)
      if (sum(X.init.sd==0)) {
        stop("at least one column in X is a constant, please remove the constant(s)")
      }
      y.init <- y
      y.init.sd <- apply(y.init,2,sd)
      y.init.mean <- mean(y.init)
      X <- scale(X,center=TRUE,scale=X.init.sd)
      y <- scale(y,center=y.init.mean,scale=y.init.sd)

      # ---- auto-dispatch ----------------------------------------------------
      # `approx = "auto"` (the default) uses the exact path when the
      # sample size doesn't exceed the would-be landmark count, and
      # switches to the Nystrom approximation otherwise (with a
      # one-line notice so the choice is visible at the call site).
      # "none" forces the exact path; "nystrom" forces the approximation.
      #
      # Note: derivative=TRUE with vcov=FALSE is *supported* under
      # Nystrom (returns point estimates without SEs), so it does
      # NOT fall through to the exact path here -- it would just
      # produce a confusing "vcovmatc not found" error there.
      if (approx == "auto") {
        if (whichkernel != "gaussian" || !is.null(eigtrunc)) {
          # Conditions truly unsupported under approx = "nystrom"
          # fall through to the exact path silently.
          approx <- "none"
        } else {
          # Validate the user-supplied nystrom_m up front so bad
          # values (NA, Inf, non-integer) produce a clean error
          # rather than being silently coerced into the dispatch
          # threshold or printed into the auto-switch message.
          effective_m <- if (is.null(nystrom_m)) {
            min(500L, n)
          } else {
            .validate_nystrom_m(nystrom_m, n)
          }
          if (n > effective_m) {
            approx <- "nystrom"
            # message() so users still see the dispatch decision under
            # print.level = 0; suppress with suppressMessages() if
            # truly silent operation is needed.
            message("krls: N = ", n, " exceeds nystrom_m = ",
                    effective_m, "; using approx = \"nystrom\". ",
                    "Set approx = \"none\" to force the exact path.")
          } else {
            approx <- "none"
          }
        }
      }

      # Re-run the derivative+vcov guard now that approx is resolved.
      # The exact path requires vcovmatc to compute derivative SEs,
      # so derivative=TRUE without vcov is invalid there. Nystrom
      # does NOT have this constraint (it returns point estimates
      # only when vcov=FALSE), so the guard only fires when we
      # actually land on the exact path.
      if (derivative == TRUE && vcov == FALSE && approx == "none") {
        stop("derivative=TRUE requires vcov=TRUE under the exact path; ",
             "set approx = \"nystrom\" for point-estimate derivatives ",
             "without SEs")
      }

      # ---- Nystrom approximation path ---------------------------------------
      # An explicit low-rank approximation. When vcov=TRUE, inference is
      # conditional on the selected landmarks and low-rank feature map.
      if (approx == "nystrom") {
        if (whichkernel != "gaussian")
          stop("approx = 'nystrom' currently supports whichkernel = 'gaussian' only")
        if (!is.null(eigtrunc))
          stop("eigtrunc is not used under approx = 'nystrom'")

        # X is already standardized here; pass the original X's centers
        # (column means of X.init) and scales (X.init.sd) so matrix-form
        # landmarks supplied in original X-scale get standardized
        # consistently.
        landmarks_resolved <- .resolve_landmarks(landmarks, landmark_method,
                                                 nystrom_m, X,
                                                 colMeans(X.init),
                                                 X.init.sd,
                                                 landmark_seed)
        noisy <- print.level > 2
        nys <- .fit_krls_nystrom(X, y, sigma, landmarks_resolved,
                                 lambda, nystrom_eps,
                                 L, U, tol, noisy,
                                 compute_vcov = vcov,
                                 lambda_method = lambda_method)

        if (print.level > 1)
          cat("Lambda that minimizes",
              if (lambda_method == "gcv") "GCV" else "Loo-Loss",
              "is:", round(nys$lambda, 5), "\n")

        yfitted_std <- nys$fitted_std
        alpha       <- nys$coeffs                       # m x 1
        lambda      <- nys$lambda
        vcovmatc    <- nys$vcov_alpha_std

        # Pointwise derivatives (continuous formula; fdskrls() handles binary).
        avgderiv <- varavgderivmat <- derivmat <- NULL
        if (derivative) {
          derivmat <- .nystrom_derivatives(X, nys$landmarks, nys$C,
                                           as.vector(alpha), sigma)
          colnames(derivmat) <- colnames(X)
          avgderiv <- matrix(colMeans(derivmat), nrow = 1)
          colnames(avgderiv) <- colnames(X)
          varavgderivmat <- .nystrom_ame_variances(X, nys$landmarks, nys$C,
                                                   vcovmatc, sigma)
          if (!is.null(varavgderivmat)) {
            colnames(varavgderivmat) <- colnames(X)
          }
          # Rescale to original units.
          derivmat <- scale(y.init.sd * derivmat,
                            center = FALSE, scale = X.init.sd)
          attr(derivmat, "scaled:scale") <- NULL
          avgderiv <- scale(as.matrix(y.init.sd * avgderiv),
                            center = FALSE, scale = X.init.sd)
          attr(avgderiv, "scaled:scale") <- NULL
          if (!is.null(varavgderivmat)) {
            varavgderivmat <- (y.init.sd / X.init.sd)^2 * varavgderivmat
            attr(varavgderivmat, "scaled:scale") <- NULL
          }
        }

        # Rescale fitted, Looe, R2 to original units. Looe matches the
        # exact path's convention: sum-of-squared standardized LOO
        # residuals times y.init.sd.
        yfitted <- as.matrix(yfitted_std * y.init.sd + y.init.mean)
        Looe    <- nys$Looe_std * y.init.sd
        R2      <- 1 - (var(y.init - yfitted) / (y.init.sd^2))

        binaryindicator <- matrix(FALSE, 1, d)
        colnames(binaryindicator) <- colnames(X)

        z <- list(
          K                  = NULL,           # not stored: memory win
          coeffs             = alpha,
          Looe               = Looe,
          fitted             = yfitted,
          X                  = X.init,
          y                  = y.init,
          sigma              = sigma,
          lambda             = lambda,
          R2                 = R2,
          derivatives        = derivmat,
          avgderivatives     = avgderiv,
          var.avgderivatives = varavgderivmat,
          vcov.c             = if (vcov) (y.init.sd^2) * vcovmatc else NULL,
          vcov.fitted        = NULL,
          binaryindicator    = binaryindicator,
          lambda_method      = lambda_method,
          # Nystrom-specific fields
          approx             = "nystrom",
          landmarks          = nys$landmarks,
          landmark_indices   = nys$landmark_indices,
          landmark_method    = nys$landmark_method,
          nystrom_m          = nys$nystrom_m,
          nystrom_eps        = nys$nystrom_eps,
          W_eigen            = nys$W_eigen,
          Dinvsqrt           = nys$Dinvsqrt,
          Sigma2             = nys$Sigma2,
          floored_count      = nys$floored_count,
          D_min_raw          = nys$D_min_raw,
          D_max_raw          = nys$D_max_raw,
          inference          = if (vcov) "conditional_nystrom" else "none"
        )
        class(z) <- "krls"

        if (derivative && binary) {
          z <- fdskrls(z)
        }

        if (print.level > 0 && derivative) {
          output <- setNames(as.vector(z$avgderivatives),
                             colnames(z$avgderivatives))
          cat("\n Average Marginal Effects:\n \n")
          print(output)
          cat("\n Quartiles of Marginal Effects:\n \n")
          print(apply(z$derivatives, 2, quantile, probs = c(.25, .5, .75)))
        }

        return(z)
      }

      # kernel matrix
      K <- NULL
      if(whichkernel=="gaussian"){ K <- gausskernel(X,sigma=sigma)}
      if(whichkernel=="linear"){K <- tcrossprod(X)}
      if(whichkernel=="poly2"){K <- (tcrossprod(X)+1)^2}
      if(whichkernel=="poly3"){K <- (tcrossprod(X)+1)^3}
      if(whichkernel=="poly4"){K <- (tcrossprod(X)+1)^4}
      if(is.null(K)){stop("No valid Kernel specified")}

      # eigenvalue decomposition
      Eigenobject <- eigen(K, symmetric = TRUE)

      # default lambda is chosen by leave-one-out optimization (golden section search)
      if (is.null(lambda)) {
       noisy <- print.level > 2
       lambda <- lambdasearch(L=L,U=U,y=y,Eigenobject=Eigenobject,tol=tol,eigtrunc=eigtrunc,noisy=noisy,lambda_method=lambda_method)

       if(print.level>1) {
         cat("Lambda that minimizes",
             if (lambda_method == "gcv") "GCV" else "Loo-Loss",
             "is:", round(lambda, 5), "\n")
       }

       } else {  # check user specified lambda
         stopifnot(is.vector(lambda),
                   length(lambda)==1,
                   is.numeric(lambda),
                   lambda>0)
      }
      # solve given LOO optimal or user specified lambda
      out <- solveforc(y = y, Eigenobject = Eigenobject, lambda = lambda, eigtrunc = eigtrunc)

      # fitted values
      yfitted <- K %*% out$coeffs

      ## var-covar matrix for c
      if (vcov == TRUE) {
        # sigma squared
        sigmasq <- as.vector((1/n) * crossprod(y - yfitted))

        if (is.null(eigtrunc)) {
          vcovmatc <- tcrossprod(multdiag(X = Eigenobject$vectors,
                                          d = sigmasq * (Eigenobject$values + lambda)^-2),
                                 Eigenobject$vectors)
        } else {
          # eigentruncation: keep only eigenvectors at least 'eigtrunc' times as large as the largest
          lastkeeper <- max(which(Eigenobject$values >= eigtrunc * Eigenobject$values[1]))
          vcovmatc <- tcrossprod(multdiag(X = Eigenobject$vectors[, 1:lastkeeper, drop = FALSE],
                                          d = sigmasq * (Eigenobject$values[1:lastkeeper] + lambda)^-2),
                                 Eigenobject$vectors[, 1:lastkeeper, drop = FALSE])
        }

        # var-covar for y hats
        vcovmatyhat <- crossprod(K, vcovmatc %*% K)

      } else {
        vcov.c      <- NULL
        vcov.fitted <- NULL
      }

      # compute derivatives
      avgderiv <- varavgderivmat <- derivmat <- NULL

 			if(derivative==TRUE){
        if(whichkernel!="gaussian"){
          stop("derivatives are only available when whichkernel='gaussian' is specified")
        }

      derivmat<-matrix(NA,n,d)
      varavgderivmat<- avgderivmat <- matrix(NA,1,d)

      rows <- cbind(rep(1:nrow(X), each = nrow(X)), 1:nrow(X))
      distances <- X[rows[,1],] - X[ rows[,2],]    # d by n*n matrix of pairwise distances
      colnames(derivmat)       <- colnames(X)
      colnames(varavgderivmat) <- colnames(X)

      for(k in 1:d){
       if(d==1){
             distk <-  matrix(distances,n,n,byrow=TRUE)
         } else {
             distk <-  matrix(distances[,k],n,n,byrow=TRUE)
        }
         L <-  distk*K
         # pointwise derivatives
         derivmat[,k] <- (-2/sigma)*L%*%out$coeff
         # variance for average derivative.
         # Identity: sum(L' V L) = (L 1)' V (L 1) where 1 is the all-ones
         # vector. Reduces this from O(n^3) to O(n^2) per predictor without
         # materializing the full L' V L product.
         r <- rowSums(L)
         varavgderivmat[1,k] <- (1/n^2) * (-2/sigma)^2 *
           as.numeric(crossprod(r, vcovmatc %*% r))
      }
       # avg derivatives
       avgderiv <- matrix(colMeans(derivmat),nrow=1)
       colnames(avgderiv) <- colnames(X)
      # get back to scale
      derivmat <- scale(y.init.sd*derivmat,center=FALSE,scale=X.init.sd)
      attr(derivmat,"scaled:scale")<- NULL
      avgderiv <- scale(as.matrix(y.init.sd*avgderiv),center=FALSE,scale=X.init.sd)
      attr(avgderiv,"scaled:scale")<- NULL
      varavgderivmat <- (y.init.sd/X.init.sd)^2*varavgderivmat
      attr(varavgderivmat,"scaled:scale")<- NULL
      }

      # get back to scale
      yfitted     <- yfitted*y.init.sd+y.init.mean
       if(vcov==TRUE){
           vcov.c      <- (y.init.sd^2)*vcovmatc
           vcov.fitted <- (y.init.sd^2)*vcovmatyhat
       } else {
         vcov.c      <- NULL
         vcov.fitted <- NULL
      }
      Looe        <- out$Le*y.init.sd
      # R square
      R2 <- 1-(var(y.init-yfitted)/(y.init.sd^2))

      # indicator for binary predictors
      binaryindicator=matrix(FALSE,1,d)
      colnames(binaryindicator) <- colnames(X)

      # return
   z <- list(K=K,
             coeffs=out$coeffs,
             Looe=Looe,
             fitted=yfitted,
             X=X.init,
             y=y.init,
             sigma=sigma,
             lambda=lambda,
             R2 = R2,
             derivatives=derivmat,
             avgderivatives=avgderiv,
             var.avgderivatives=varavgderivmat,
             vcov.c=vcov.c,
             vcov.fitted=vcov.fitted,
             binaryindicator=binaryindicator,
             lambda_method=lambda_method
            )
  class(z) <- "krls"

  # add first differences if requested
    if(derivative==TRUE && binary==TRUE){
      z <- fdskrls(z)
    }

  # printing
      if(print.level>0 && derivative==TRUE){
      output <- setNames(as.vector(z$avgderivatives), colnames(z$avgderivatives))
      cat("\n Average Marginal Effects:\n \n")
      print(output)

      cat("\n Quartiles of Marginal Effects:\n \n")
      print(apply(z$derivatives,2,quantile,probs=c(.25,.5,.75)))
      }

  return(z)

}
