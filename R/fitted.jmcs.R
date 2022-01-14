##' @title Fitted values for joint models
##' @name fitted
##' @aliases fitted.jmcs
##' @description Extract fitted values for joint models.
##' @param object an object inheriting from class \code{jmcs}.
##' @param type for which type of fitted values to calculate.
##' @param process for which sub-model to calculate the fitted values.
##' @param ... further arguments passed to or from other methods.
##' @return a numeric vector of fitted values.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @examples
##' \donttest{
##' fit <- jmcs(ydata = ydata, cdata = cdata, 
##'             long.formula = response ~ time + gender + x1 + race, 
##'             surv.formula = Surv(surv, failure_type) ~ x1 + gender + x2 + race, 
##'             random =  ~ time| ID)
##'
##' # fitted for the longitudinal process
##' head(cbind(
##'   "Marg" = fitted(fit, type = "Marginal", process = "Longitudinal"), 
##'   "Subj" = fitted(fit, type = "Subject", process = "Longitudinal")
##' ))
##' # fitted for the levent process - marginal survival function
##' head(fitted(fit, type = "Marginal", process = "Event"))
##' }
##' 
##' @export
##' 

fitted.jmcs <- function(object, type = c("Marginal", "Subject"), 
                        process = c("Longitudinal", "Event"), ...) {
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  
  if (type == "Marginal" & process == "Longitudinal") {
    fitted <- object$fitted$fittedmar 
  } else if (type == "Subject" & process == "Longitudinal") {
    fitted <- object$fitted$fitted 
  } else if (type == "Marginal" & process == "Event") {
    fitted <- object$fittedSurv
  } else {
    stop("Please choose one of the following options: Marginal, Subject.")
  }
  fitted
}