##' @title Extract log-likelihoods from fitted joint models
##' @name logLik
##' @description Extract the stored, numerically evaluated marginal
##'   log-likelihood from a fitted single-failure or competing-risks model.
##' @param object A fitted model of class \code{jmcs}, \code{JMMLSM},
##'   or \code{mvjmcs}.
##' @param ... further arguments passed to or from other methods.
##' @details The likelihood is returned without rounding or recomputation.
##'   The parameter-count convention includes regression coefficients,
##'   association parameters, variance parameters, and one baseline-hazard
##'   increment per row of \code{H01} (and \code{H02} for competing risks).
##'   Random-effects covariance matrices contribute their unique elements;
##'   subject-specific random-effect predictions are not counted.
##'   The \code{nobs} attribute is the number of independent subjects,
##'   represented by survival-data rows for \code{jmcs} and \code{JMMLSM}.
##'   For \code{mvjmcs}, it uses the stored fitted-subject count, which
##'   accounts for landmark filtering rather than counting the raw data.
##'
##'   This is a full parameter-count convention for the estimated
##'   baseline-hazard increments, not a claim about their effective degrees
##'   of freedom. Consequently, AIC and BIC computed from this object use
##'   this explicit counting convention. Compare models fitted to the same
##'   data with compatible likelihood definitions and integration settings.
##'
##'   For \code{mvjmcs}, the value is the stored Laplace approximation to
##'   the marginal log-likelihood of the fitted working model. It is not
##'   exact integration, and the fitting algorithm is not guaranteed to
##'   maximize this reported approximation exactly. Landmark results refer
##'   to the retained fitting sample and do not introduce a separate
##'   selection-normalized likelihood. Older fitted objects without a
##'   stored likelihood and fitted-subject count must be refitted.
##' @return A scalar of class \code{logLik}, with attributes \code{df}
##'   (parameter count) and \code{nobs} (number of subjects). For
##'   \code{mvjmcs}, the \code{approximation} attribute is \code{"Laplace"}.
##' @seealso \code{\link{jmcs}}, \code{\link{JMMLSM}},
##'   \code{\link{mvjmcs}}, \code{\link[stats]{AIC}}, \code{\link[stats]{BIC}}
##' @exportS3Method stats::logLik
logLik.jmcs <- function(object, ...) {

  if (!inherits(object, "jmcs")) stop("Use only with 'jmcs' objects.")
  ll <- object$loglike
  if (!is.numeric(ll) || length(ll) != 1L || !is.finite(ll))
    stop("The fitted object has no finite scalar 'loglike' value.")
  if (!is.null(object$convergence) &&
      (length(object$convergence) != 1L || is.na(object$convergence) ||
       object$convergence != 1))
    stop("The fitted model did not converge.")
  if (!is.logical(object$CompetingRisk) || length(object$CompetingRisk) != 1L ||
      is.na(object$CompetingRisk)) stop("Invalid 'CompetingRisk' flag.")

  parameters <- c("beta", "gamma1", "nu1", "sigma")
  if (object$CompetingRisk) parameters <- c(parameters, "gamma2", "nu2")
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
