# Internal: Laplace marginal log-likelihood for the working model used
# by getbSig/getbSigSF. Call only after the final E-step, at final estimates.
# Those objectives include the normal random-effects constant, but omit
# the longitudinal N_i * log(2*pi)/2 constant, which is restored below.
getmvLoglike <- function(data, posterior) {
  n <- length(data$Y)
  if (!n || length(posterior) != n)
    stop("Posterior results must match the subjects in the fitting data.")
  contributions <- rep(NA_real_, n)
  failed <- integer(0)
  for (i in seq_len(n)) {
    fit <- posterior[[i]]
    q <- length(fit$mode)
    U <- fit$ccov
    status <- fit$optim.convergence
    ok <- length(status) == 1L && !is.na(status) && status == 0L &&
      is.numeric(fit$objective) && length(fit$objective) == 1L &&
      is.finite(fit$objective) && q > 0L &&
      is.matrix(U) && identical(dim(U), c(q, q)) &&
      all(is.finite(U)) && all(diag(U) > 0)
    if (!ok) {
      failed <- c(failed, i)
      next
    }
    n_long <- sum(vapply(data$Y[[i]], length, integer(1)))
    # crossprod(U) = inverse Hessian, hence sum(log(diag(U)))
    # equals -log(det(Hessian))/2.
    contributions[i] <- -as.numeric(fit$objective) +
      (q - n_long) * log(2 * pi) / 2 + sum(log(diag(U)))
    if (!is.finite(contributions[i])) failed <- c(failed, i)
  }
  if (length(failed))
    warning("Laplace log-likelihood unavailable: invalid final posterior ",
            "optimization for subject indices ", paste(failed, collapse = ", "),
            ".", call. = FALSE)
  list(value = if (length(failed)) NA_real_ else sum(contributions),
       contributions = contributions, nobs = n,
       method = "Laplace", failed = failed)
}
