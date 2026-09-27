##' @title Extract parameter estimates from fitted joint models
##' @name coef
##' @description Extract unrounded regression, association, and
##'   variance parameters from fitted single-failure or competing-risks
##'   joint models of class \code{jmcs}, \code{JMMLSM}, or \code{mvjmcs}.
##' @param object A fitted model of class \code{jmcs},
##'   \code{JMMLSM}, or \code{mvjmcs}.
##' @param ... Further arguments passed to or from other methods.
##' @details
##' Parameters are returned in the order used by the corresponding
##' \code{vcov()} method:
##' \describe{
##'   \item{\code{jmcs}}{
##'     Longitudinal regression coefficients, event regression
##'     coefficients, association parameters, residual variance,
##'     and random-effects covariance parameters.
##'   }
##'   \item{\code{JMMLSM}}{
##'     Longitudinal mean and variance regression coefficients,
##'     event regression coefficients, mean association parameters,
##'     variance association parameters, and random-effects
##'     covariance parameters.
##'   }
##'   \item{\code{mvjmcs}}{
##'     Longitudinal regression coefficients across biomarkers,
##'     residual variances for each biomarker, event regression
##'     coefficients, association parameters, and random-effects
##'     covariance parameters. The association structures
##'     \code{"sre"}, \code{"present"}, and \code{"presentlp"}
##'     are supported.
##'   }
##' }
##' Within each event or association parameter group, parameters
##' for event 1 precede those for event 2. Event-2 parameters are
##' included only for competing-risks models.
##'
##' Random-effects variances precede covariances. Covariances are
##' ordered by increasing distance from the diagonal.
##' Baseline hazard increments and subject-specific random-effect
##' predictions are not included.
##' @return A named numeric vector with names and ordering matching
##'   the rows and columns of \code{vcov(object)}.
##'   Estimates are returned without rounding or transformation:
##'   variance components remain variances, and variance-model
##'   regression coefficients are not exponentiated.
##' @seealso \code{\link{jmcs}}, \code{\link{JMMLSM}},
##'   \code{\link{mvjmcs}}, \code{\link[stats]{vcov}}
##' @export

coef.jmcs <- function(object, ...) {
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' objects.")
  if (!is.logical(object$CompetingRisk) ||
      length(object$CompetingRisk) != 1L || is.na(object$CompetingRisk))
    stop("The fitted object must contain a valid 'CompetingRisk' flag.")
  
  required <- c("beta", "gamma1", "nu1", "sigma", "Sig")
  if (object$CompetingRisk)
    required <- c(required, "gamma2", "nu2")
  valid <- vapply(required, function(nm) {
    x <- object[[nm]]
    is.numeric(x) && length(x) > 0L && all(is.finite(x))
  }, logical(1))
  if (!all(valid))
    stop("Missing or invalid parameter estimates: ",
         paste(required[!valid], collapse = ", "), ".")
  
  S <- as.matrix(object$Sig)
  q <- nrow(S)
  if (q < 1L || ncol(S) != q || length(object$nu1) != q ||
      length(object$sigma) != 1L ||
      (object$CompetingRisk && length(object$nu2) != q))
    stop("Inconsistent parameter dimensions in the fitted object.")
  
  named_block <- function(x, prefix, label) {
    nm <- names(x)
    if (is.null(nm) || length(nm) != length(x) ||
        anyNA(nm) || any(!nzchar(nm)))
      stop("Missing parameter names in '", label, "'.")
    stats::setNames(as.numeric(x), paste0(prefix, nm))
  }
  
  # Match the existing jmcs association labels.
  re_names <- c("(Intercept)", all.vars(object$random)[seq_len(q - 1L)])
  if (length(re_names) != q || anyNA(re_names))
    stop("Cannot identify random-effect names from 'object$random'.")
  
  ans <- c(named_block(object$beta, "Y.", "beta"),
           named_block(object$gamma1, "T.", "gamma1"))
  if (object$CompetingRisk)
    ans <- c(ans, named_block(object$gamma2, "T.", "gamma2"))
  ans <- c(ans, stats::setNames(as.numeric(object$nu1),
                                paste0("T.asso:", re_names, "_1")))
  if (object$CompetingRisk)
    ans <- c(ans, stats::setNames(as.numeric(object$nu2),
                                  paste0("T.asso:", re_names, "_2")))
  ans <- c(ans, stats::setNames(as.numeric(object$sigma), "Y.sigma^2"))
  
  # getCov.cpp/getCovSF.cpp: diagonals, then successive upper diagonals.
  covariance <- stats::setNames(diag(S), paste0("Sig", seq_len(q), seq_len(q)))
  for (lag in seq_len(q - 1L)) {
    i <- seq_len(q - lag)
    j <- i + lag
    covariance <- c(covariance,
                    stats::setNames(S[cbind(i, j)], paste0("Sig", i, j)))
  }
  c(ans, covariance)
}
