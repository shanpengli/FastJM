##' @title Residuals for joint models
##' @name residuals
##' @aliases residuals.jmcs
##' @description Extract residuals for joint models.
##' @param object an object inheriting from class \code{jmcs}.
##' @param type what type of residuals to calculate. 
##' @param ... further arguments passed to or from other methods.
##' @return a vector of residuals of the longitudinal sub-model.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}}
##' @examples 
##' \donttest{
##' # a joint model fit
##' fit <- jmcs(ydata = ydata, cdata = cdata, 
##' long.formula = response ~ time + x1, 
##' surv.formula = Surv(surv, failure_type) ~ x1 + x2, 
##' random =  ~ time| ID)
##' 
##' # residuals of the longitudinal sub-model
##' head(cbind(
##'   "Marg" = residuals(fit, type = "Marginal"), 
##'   "Subj" = residuals(fit, type = "Subject")
##' ))
##' }
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