##' @title Variance-covariance matrix of the estimated parameters for joint models
##' @name vcov
##' @aliases vcov.jmcs
##' @description Extract variance-covariance matrix for joint models.
##' @param object an object inheriting from class \code{jmcs}.
##' @param ... further arguments passed to or from other methods.
##' @return a matrix of variance covariance of all parameter estimates.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}}
##' @export
##' 

vcov.jmcs <- function(object, ...) {
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  
  vcov <- as.data.frame(object$vcov)
  long <- all.vars(object$LongitudinalSubmodel)[-1]
  survival <- all.vars(object$SurvivalSubmodel)[-(1:2)]
  survival1 <- paste0("T.", survival, "_1")
  long <- paste0("Y.", c("(Intercept)", long))
  p1a <- length(object$nu1)
  nu1 <- rep("nu1", p1a)
  for (i in 1:p1a) nu1[i] <- paste0("T.", nu1[i], "_", i)
  sig <- "Sig"
  if (p1a == 1) {
    sig <- paste0(sig, "11") 
  } else {
    sig <- rep(sig, p1a*(p1a+1)/2)
    for (i in 1:p1a) sig[i] <- paste0(sig[i], i, i)
    if (p1a == 2) sig[p1a+1] <- paste0(sig[p1a+1], "12")
    if (p1a == 3) {
      sig[p1a+1] <- paste0(sig[p1a+1], "12")
      sig[p1a+2] <- paste0(sig[p1a+2], "23")
      sig[p1a+3] <- paste0(sig[p1a+3], "13")
    }
  }
  
  if (object$CompetingRisk) {
    survival2 <- paste0("T.", survival, "_2")
    nu2 <- rep("nu2", p1a)
    for (i in 1:p1a) nu2[i] <- paste0("T.", nu2[i], "_", i)
    colnames(vcov) <- c(long, survival1, survival2, nu1, nu2, "Y.sigma^2", sig)
    rownames(vcov) <- c(long, survival1, survival2, nu1, nu2, "Y.sigma^2", sig)
  } else {
    colnames(vcov) <- c(long, survival1, nu1, "Y.sigma^2", sig)
    rownames(vcov) <- c(long, survival1, nu1, "Y.sigma^2", sig)
  }
  
  vcov
}

