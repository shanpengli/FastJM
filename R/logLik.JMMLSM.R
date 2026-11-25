##' @rdname logLik
##' @exportS3Method stats::logLik
logLik.JMMLSM <- function(object, ...) {

  if (!inherits(object, "JMMLSM")) stop("Use only with 'JMMLSM' objects.")
  ll <- object$loglike
  if (!is.numeric(ll) || length(ll) != 1L || !is.finite(ll))
    stop("The fitted object has no finite scalar 'loglike' value.")
  if (!is.null(object$convergence) &&
      (length(object$convergence) != 1L || is.na(object$convergence) ||
       object$convergence != 1))
    stop("The fitted model did not converge.")
  if (!is.logical(object$CompetingRisk) || length(object$CompetingRisk) != 1L ||
      is.na(object$CompetingRisk)) stop("Invalid 'CompetingRisk' flag.")

  parameters <- c("beta", "tau", "gamma1", "alpha1", "vee1")
  if (object$CompetingRisk) parameters <- c(parameters, "gamma2", "alpha2", "vee2")
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

  cd <- object$cdata
  if (!(is.data.frame(cd) || is.matrix(cd)) || nrow(cd) < 1L)
    stop("Missing survival data needed to determine the number of subjects.")
  vars <- all.vars(object$random)
  if (!length(vars)) stop("Missing subject identifier in the random-effects formula.")
  id <- vars[length(vars)]
  if (!id %in% colnames(cd) || anyNA(cd[, id]) || anyDuplicated(cd[, id]))
    stop("Survival data must contain one row per subject with a valid identifier.")
  structure(as.numeric(ll), df = as.numeric(df), nobs = nrow(cd), class = "logLik")
}
