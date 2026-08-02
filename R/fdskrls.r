fdskrls <-
  function(object,...){
    d <- ncol(object$X)
    n <- nrow(object$X)
    lengthunique    <- function(x){length(unique(x))}

    # SEs are only available when vcov was computed. If unavailable
    # (for example vcov=FALSE), still return FD point estimates.
    has_se <- !is.null(object$vcov.c)

    fdderivatives           =object$derivatives
    fdavgderivatives        =object$avgderivatives
    fdvar.var.avgderivatives=object$var.avgderivatives

    # Positions of binary variables.
    #
    # Two different spaces are involved. The counterfactual data X1/X0 below is
    # built from `object$X`, the RAW input. The matrices written at the end --
    # derivatives, avgderivatives, binaryindicator -- live in the PROCESSED
    # (kernel) design. When cat_columns one-hot expands a factor the two stop
    # lining up, so a single index vector cannot serve both: raw positions
    # written into processed matrices land on the wrong column.
    #
    # `prep$col_origin` gives, for each processed column, the raw column it came
    # from, so it resolves the two spaces without relying on column names (which
    # are absent when X is an unnamed matrix).
    binaryindicator <- which(apply(object$X, 2, lengthunique) == 2)
    col.origin <- object$prep$col_origin
    if (!is.null(col.origin) && is.matrix(object$derivatives) &&
        length(col.origin) == ncol(object$derivatives)) {
      # A raw binary column that was itself expanded (e.g. a two-level factor
      # named in cat_columns) has no single processed counterpart; leave it to
      # the kernel derivative rather than writing to an arbitrary column.
      proc.index <- vapply(binaryindicator, function(j) {
        pos <- which(col.origin == j)
        if (length(pos) == 1L) pos else NA_integer_
      }, integer(1))
      keep <- !is.na(proc.index)
      binaryindicator <- binaryindicator[keep]
      proc.index      <- proc.index[keep]
    } else {
      proc.index <- binaryindicator
    }
      if(length(binaryindicator)==0){
        # no binary vars in X; return derivs as is
      } else {
        # compute marginal differences from min to max
        est  <- se <- matrix(NA,nrow=1,ncol=length(binaryindicator))
        diffsstore <- matrix(NA,nrow=n,ncol=length(binaryindicator))
        for(i in 1:length(binaryindicator)){
          X1 <- X0 <- object$X
          # test data with D=Max
          X1[,binaryindicator[i]] <- max(X1[,binaryindicator[i]])
          # test data with D=Min
          X0[,binaryindicator[i]] <- min(X0[,binaryindicator[i]])
          Xall      <- rbind(X1,X0)
          # contrast vector
          h         <- matrix(rep(c(1/n,-(1/n)),each=n),ncol=1)
          # fitted values
          pout      <- predict(object,newdata=Xall,se.fit = has_se)
          # store FD estimates
          est[1,i] <- t(h)%*%pout$fit
          if (has_se) {
            # SE (multiply by sqrt2 to correct for using data twice )
            se[1,i] <- as.vector(sqrt(t(h)%*%pout$vcov.fit%*%h))*sqrt(2)
          }
          # all
          diffs <- pout$fit[1:n]-pout$fit[(n+1):(2*n)]
          diffsstore[,i] <- diffs
        }
        # sub in first differences, at the PROCESSED column positions
        object$derivatives[,proc.index] <- diffsstore
        object$avgderivatives[,proc.index] <- est
        if (has_se && !is.null(object$var.avgderivatives)) {
          object$var.avgderivatives[,proc.index] <- se^2
        }
        object$binaryindicator[,proc.index] <- TRUE
      }

   return(invisible(object))

}
