##' Control options for \code{JMMLSM()}
##'
##' Specifies numerical integration, convergence, initialization, and optimizer
##' settings for \code{\link{JMMLSM}()}.
##'
##' @title Control Options for JMMLSM
##'
##' @param maxiter Maximum number of EM iterations. The default is \code{1000}.
##' @param tol Convergence tolerance. The default is \code{1e-4}.
##' @param quadpoint Number of Gauss--Hermite quadrature points. The default is
##' \code{15}.
##' @param verbose Logical value indicating whether to print iteration details.
##' The default is \code{FALSE}.
##' @param initial.para Optional list of user-supplied initial parameter values.
##' The default is \code{NULL}.
##' @param method Numerical integration method used in the E-step. Available
##' options are \code{"adaptive"} and \code{"standard"}. The default is
##' \code{"adaptive"}.
##' @param opt Optimization method used to fit the initial linear mixed-effects
##' model. Available options are \code{"nlminb"} and \code{"optim"}. The default
##' is \code{"nlminb"}.
##'
##' @return A list of control parameters used by \code{\link{JMMLSM}()}.
##' @export
JMMLSM_control <- function(maxiter = 1000,
                           tol = 1e-4,
                           quadpoint = NULL,
                           verbose = FALSE,
                           initial.para = NULL,
                           method = c("adaptive", "standard"),
                           opt = c("nlminb", "optim")) {
  
  method <- match.arg(method)
  opt <- match.arg(opt)
  
  if (!is.numeric(maxiter) || length(maxiter) != 1 || maxiter <= 0) {
    stop("'maxiter' must be a positive number.")
  }
  
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("'tol' must be a positive number.")
  }
  
  if (!is.null(quadpoint) &&
      (!is.numeric(quadpoint) || length(quadpoint) != 1 || quadpoint <= 0)) {
    stop("'quadpoint' must be a positive number or NULL.")
  }
  
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be TRUE or FALSE.")
  }
  
  list(
    maxiter = maxiter,
    tol = tol,
    quadpoint = quadpoint,
    verbose = verbose,
    initial.para = initial.para,
    method = method,
    opt = opt
  )
}