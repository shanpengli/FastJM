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
  
  getdum <- getdummy(long.formula = object$LongitudinalSubmodel, 
                     surv.formula = object$SurvivalSubmodel, 
                     random = object$random, ydata = object$ydata, cdata = object$cdata)
  
  surv.formula <- getdum$surv.formula
  random <- all.vars(object$random)
  vcov <- as.data.frame(object$vcov)
  long <- names(object$beta)
  survival <- names(object$gamma1)
  survival1 <- paste0("T.", survival)
  long <- paste0("Y.", long)
  p1a <- length(object$nu1)
  nu1 <- rep("T.asso:", p1a)
  sig <- "Sig"
  if (p1a == 1) {
    sig <- paste0(sig, "11")
    nu1 <- paste0(nu1, "(Intercept)_1")
  } else {
    sig <- rep(sig, p1a*(p1a+1)/2)
    for (i in 1:p1a) sig[i] <- paste0(sig[i], i, i)
    if (p1a == 2) {
      sig[p1a+1] <- paste0(sig[p1a+1], "12")
      nu1[1] <- paste0(nu1[1], "(Intercept)_1")
      nu1[2] <- paste0(nu1[2], random[1], "_1")
    }
    if (p1a == 3) {
      sig[p1a+1] <- paste0(sig[p1a+1], "12")
      sig[p1a+2] <- paste0(sig[p1a+2], "23")
      sig[p1a+3] <- paste0(sig[p1a+3], "13")
      nu1[1] <- paste0(nu1[1], "(Intercept)_1")
      nu1[2] <- paste0(nu1[2], random[1], "_1")
      nu1[3] <- paste0(nu1[3], random[2], "_1")
    }
  }
  
  if (object$CompetingRisk) {
    survival <- names(object$gamma2)
    survival2 <- paste0("T.", survival)
    nu2 <- rep("T.asso:", p1a)
    if (p1a == 1) nu2 <- paste0(nu2, "(Intercept)_2")
    if (p1a == 2) {
      nu2[1] <- paste0(nu2[1], "(Intercept)_2")
      nu2[2] <- paste0(nu2[2], random[1], "_2")
    }
    if (p1a == 3) {
      nu2[1] <- paste0(nu2[1], "(Intercept)_2")
      nu2[2] <- paste0(nu2[2], random[1], "_2")
      nu2[3] <- paste0(nu2[3], random[2], "_2")
    }

    colnames(vcov) <- c(long, survival1, survival2, nu1, nu2, "Y.sigma^2", sig)
    rownames(vcov) <- c(long, survival1, survival2, nu1, nu2, "Y.sigma^2", sig)
  } else {
    colnames(vcov) <- c(long, survival1, nu1, "Y.sigma^2", sig)
    rownames(vcov) <- c(long, survival1, nu1, "Y.sigma^2", sig)
  }
  
  vcov
}

