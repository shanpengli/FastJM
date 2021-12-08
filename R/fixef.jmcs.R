##' @title Estimated coefficients estimates for joint models
##' @name fixef
##' @aliases fixef.jmcs
##' @description Extracts the fixed effects for a fitted joint model.
##' @param object an object inheriting from class \code{jmcs}.
##' @param process for which sub-model to extract the estimated coefficients.
##' @param ... further arguments passed to or from other methods.
##' @export
##' 

fixef.jmcs <- function(object, process = c("Longitudinal", "Event"), ...) {
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  
  if (process == "Longitudinal") {
    vars <- object$beta
  } else if (process == "Event") {
    if (!is.null(object$gamma2)) {
      vars <- vector("list", 2)
      vars[[1]] <- object$gamma1
      vars[[2]] <- object$gamma2
      names(vars) <- c("Risk1", "Risk2")
    } else {
      vars <- vector("list", 1)
      vars[[1]] <- object$gamma1
      names(vars) <- c("Risk1")
    }
  } else {
    stop("Please choose one the following arguments: Longitudinal, Event.")
  }
  vars

}