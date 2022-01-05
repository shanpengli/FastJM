##' @title Estimated coefficients estimates for joint models
##' @name fixef
##' @description Extracts the fixed effects for a fitted joint model.
##' @param object an object inheriting from class \code{jmcs}.
##' @param process for which sub-model to extract the estimated coefficients.
##' @param ... further arguments passed to or from other methods.
##' @return A numeric vector or a list of the estimated parameters for the fitted model.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @examples 
##' \donttest{
##' # a joint model fit
##' fit <- jmcs(ydata = ydata, cdata = cdata, 
##' long.formula = response ~ time + x1, 
##' surv.formula = Surv(surv, failure_type) ~ x1 + x2, 
##' random =  ~ time| ID)
##' 
##' # fixed effects for the longitudinal process
##' fixef(fit, process = "Longitudinal")
##' # fixed effects for the event process
##' fixef(fit, process = "Event")
##' }
##' @export
##' 

fixef <- function(object, process = c("Longitudinal", "Event"), ...) {
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