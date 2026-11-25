##' @rdname logLik
##' @exportS3Method stats::logLik
logLik.mvjmcs <- function(object, ...) {

  if (!inherits(object, "mvjmcs")) stop("Use only with 'mvjmcs' objects.")
  ll <- object$loglike
  if (!is.numeric(ll) || length(ll) != 1L || !is.finite(ll))
    stop("No finite stored log-likelihood. Refit with the updated mvjmcs() ",
         "and check model and posterior-optimization convergence.", call. = FALSE)
  if (!is.null(object$convergence) &&
      (length(object$convergence) != 1L || is.na(object$convergence) ||
       object$convergence != 1))
    stop("The fitted model did not converge.")
  if (!is.logical(object$CompetingRisk) || length(object$CompetingRisk) != 1L ||
      is.na(object$CompetingRisk)) stop("Invalid 'CompetingRisk' flag.")

  parameters <- c("beta", "sigma", "gamma1", "alpha1")
  if (object$CompetingRisk) parameters <- c(parameters, "gamma2", "alpha2")
  valid <- vapply(parameters, function(nm) {
    x <- object[[nm]]
    is.numeric(x) && length(x) > 0L && all(is.finite(x))
  }, logical(1))
  if (!all(valid)) stop("Missing or invalid model parameters: ",
                       paste(parameters[!valid], collapse = ", "), ".")
  S <- object$Sig
  if (!is.matrix(S) || !is.numeric(S) || nrow(S) < 1L ||
      ncol(S) != nrow(S) || any(!is.finite(S)))
    stop("Missing or invalid random-effects covariance matrix 'Sig'.")
  q <- nrow(S)
  df <- sum(vapply(parameters, function(nm) length(object[[nm]]), integer(1))) +
    q * (q + 1) / 2

  hazards <- if (object$CompetingRisk) c("H01", "H02") else "H01"
  for (nm in hazards) {
    H <- object[[nm]]
    if (!(is.matrix(H) || is.data.frame(H)) || ncol(H) != 3L ||
        !is.numeric(H[, 1L]) || any(!is.finite(H[, 1L])) ||
        anyDuplicated(H[, 1L]))
      stop("Invalid baseline-hazard table '", nm, "'.")
    df <- df + nrow(H)
  }

  # Use the fitted sample size, not raw cdata, which may include subjects
  # excluded during landmark preprocessing.
  n <- object$nobs
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n != floor(n))
    stop("Missing or invalid fitted-subject count. Refit with the updated mvjmcs().")
  if (!identical(object$loglike.method, "Laplace"))
    stop("Unrecognized or missing likelihood approximation. Refit with the updated mvjmcs().")
  structure(as.numeric(ll), df = as.numeric(df), nobs = as.integer(n),
            approximation = "Laplace", class = "logLik")
}
