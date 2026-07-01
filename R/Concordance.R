##' @title Concordance for joint models
##' @name Concordance
##' @aliases Concordance
##'
##' @param seed A numeric value used to set the random seed for cross-validation.
##' Default is \code{100}.
##' @param object A fitted object of class \code{jmcs} or \code{JMMLSM}.
##' @param n.cv Number of cross-validation folds. Default is \code{3}.
##' @param maxiter Maximum number of EM iterations allowed when refitting the
##' model within each cross-validation fold. Default is \code{10000}.
##' @param initial.para Logical value indicating whether initial parameter values
##' are used when refitting the model. Default is \code{TRUE}.
##' @param ... Further arguments passed to model-specific methods.
##'
##' @return
##' A list containing:
##' \describe{
##'   \item{\code{n.cv}}{The number of cross-validation folds.}
##'   \item{\code{Concordance.cv}}{A list of fold-specific concordance estimates.}
##'   \item{\code{PI.cv}}{A list of fold-specific Cindex value for the estimated time-independent prognostic indexes across subjects.}
##'   \item{\code{CompetingRisk}}{Logical value indicating whether the fitted model accounts for competing risks.}
##'   \item{\code{seed}}{The random seed used for cross-validation.}
##' }
##'
##' @seealso \code{\link{jmcs}}, \code{\link{JMMLSM}}
##' @export
Concordance <- function(seed = 100, object, n.cv = 3, maxiter = 10000,
                        initial.para = TRUE, ...) {
  
  if (!any(class(object) %in% c("jmcs", "JMMLSM"))) {
    stop("Use only with 'jmcs' or 'JMMLSM' objects.\n")
  }
  
  if (inherits(object, "jmcs")) {
    res <- Concordance.jmcs(seed = seed,
                            object = object,
                            n.cv = n.cv,
                            maxiter = maxiter,
                            initial.para = initial.para,
                            ...)
  } else {
    res <- Concordance.JMMLSM(seed = seed,
                              object = object,
                              n.cv = n.cv,
                              maxiter = maxiter,
                              initial.para = initial.para,
                              ...)
  }
  
  return(res)
}