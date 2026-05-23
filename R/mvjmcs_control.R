##' Control options for \code{mvjmcs()}
##'
##' Specifies convergence, initialization, optimizer, verbosity, and parallel
##' computing settings for \code{\link{mvjmcs}()}.
##'
##' @title Control Options for mvjmcs
##'
##' @param maxiter Maximum number of EM iterations. The default is \code{10000}.
##' @param opt Optimization method used to fit the initial mixed-effects models.
##' Available options are \code{"nlminb"} and \code{"optim"}. The default is
##' \code{"nlminb"}.
##' @param tol Convergence tolerance. The default is \code{0.005}.
##' @param verbose Logical value indicating whether to print iteration details.
##' The default is \code{FALSE}.
##' @param initial.para Optional list of user-supplied initial parameter values
##' for the EM algorithm. The default is \code{NULL}.
##' @param cpu.cores Number of CPU cores used for parallel computation. The
##' default is \code{NULL}, in which case the function uses its default
##' computation strategy.
##'
##' @return A list of control parameters used by \code{\link{mvjmcs}()}.
##'
##' @export
mvjmcs_control <- function(maxiter = 10000,
                           opt = c("nlminb", "optim"),
                           tol = 0.005,
                           verbose = FALSE,
                           initial.para = NULL,
                           cpu.cores = NULL) {
  
  opt <- match.arg(opt)
  
  if (!is.numeric(maxiter) || length(maxiter) != 1 || maxiter <= 0) {
    stop("'maxiter' must be a positive number.")
  }
  
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("'tol' must be a positive number.")
  }
  
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be TRUE or FALSE.")
  }
  
  if (!is.null(cpu.cores) &&
      (!is.numeric(cpu.cores) || length(cpu.cores) != 1 ||
       cpu.cores <= 0 || cpu.cores != floor(cpu.cores))) {
    stop("'cpu.cores' must be a positive integer or NULL.")
  }
  
  list(
    maxiter = maxiter,
    opt = opt,
    tol = tol,
    verbose = verbose,
    initial.para = initial.para,
    cpu.cores = cpu.cores
  )
}