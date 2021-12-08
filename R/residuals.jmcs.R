##' @title Residuals for joint models
##' @name residuals
##' @aliases residuals.jmcs
##' @description Extract residuals for joint models.
##' @param object an object inheriting from class \code{jmcs}.
##' @param type what type of residuals to calculate. 
##' @param ... further arguments passed to or from other methods.
##' @export
##' 

residuals.jmcs <- function(object, type = c("Marginal", "Subject"), ...) {
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  
  if (type == "Marginal") {
    resid <- object$fitted$residmar 
  } else if (type == "Subject") {
    resid <- object$fitted$resid 
  } else {
    stop("Please choose one of the following options: Marginal, Subject.")
  }
  resid
}