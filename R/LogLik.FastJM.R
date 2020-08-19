logLik.FastJM <-
  function (object, ...) {
    if (!inherits(object, "FastJM"))
      stop("Use only with 'FastJM' objects.\n")
    out <- object$loglike
    attr(out, "nobs") <- object$k
    class(out) <- "logLik"
    out
  }
