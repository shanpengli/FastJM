##' Control options for \code{jmcs()}
##'
##' Specifies numerical integration, convergence, initialization, and optimizer
##' settings for \code{\link{jmcs}()}.
##'
##' @title Control Options for jmcs
##'
##' @param quadpoint Number of Gauss--Hermite quadrature points used for numerical
##' integration. The default is \code{NULL}.
##' @param maxiter Maximum number of EM iterations. The default is \code{10000}.
##' @param verbose Logical value indicating whether to print iteration details.
##' The default is \code{FALSE}.
##' @param initial.para Optional list of user-supplied initial parameter values
##' for the EM algorithm. The default is \code{NULL}.
##' @param tol Convergence tolerance. The default is \code{1e-4}.
##' @param method Numerical integration method used in the E-step. Available
##' options are \code{"pseudo-adaptive"} and \code{"standard"}. The default is
##' \code{"pseudo-adaptive"}.
##' @param opt Optimization method used to fit the initial linear mixed-effects
##' model. Available options are \code{"nlminb"} and \code{"optim"}. The default
##' is \code{"nlminb"}.
##'
##' @return A list of control parameters used by \code{\link{jmcs}()}.
##'
##' @export
jmcs_control <- function(quadpoint = NULL,
                         maxiter = 10000,
                         verbose = FALSE,
                         initial.para = NULL,
                         tol = 1e-4,
                         method = c("pseudo-adaptive", "standard"),
                         opt = c("nlminb", "optim")) {
  
  method <- match.arg(method)
  opt <- match.arg(opt)
  
  if (!is.null(quadpoint) &&
      (!is.numeric(quadpoint) || length(quadpoint) != 1 || quadpoint <= 0)) {
    stop("'quadpoint' must be a positive number or NULL.")
  }
  
  if (!is.numeric(maxiter) || length(maxiter) != 1 || maxiter <= 0) {
    stop("'maxiter' must be a positive number.")
  }
  
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("'tol' must be a positive number.")
  }
  
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be TRUE or FALSE.")
  }
  
  list(
    quadpoint = quadpoint,
    maxiter = maxiter,
    verbose = verbose,
    initial.para = initial.para,
    tol = tol,
    method = method,
    opt = opt
  )
}

