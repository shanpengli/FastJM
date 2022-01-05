##' @title Random effects estimates for joint models
##' @name ranef
##' @description Extracts the posterior mean of the random effects for a fitted joint model.
##' @param object an object inheriting from class \code{jmcs}.
##' @param ... further arguments passed to or from other methods.
##' @return a matrix of random effects estimates.
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
##' # extract random effects estimates
##' head(ranef(fit))
##' }
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