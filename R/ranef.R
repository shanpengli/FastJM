##' @title Random effects estimates for joint models
##' @name ranef
##' @description Extracts the posterior mean of the random effects for a fitted joint model.
##' @param object an object inheriting from class \code{jmcs}.
##' @param ... further arguments passed to or from other methods.
##' @export
##' 

ranef <- function(object, ...) {
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  
  rand <- all.vars(object$random)
  dat <- data.frame(t(object$FUNB))
  if (length(rand) == 1) {
    colnames(dat) <- c("(Intercept)") 
  } else {
    colnames(dat) <- c("(Intercept)", rand[-length(rand)]) 
  }
  dat
}