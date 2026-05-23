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
  
  vcov <- as.data.frame(object$vcov)
  
  random <- all.vars(object$random)
  p1a <- length(object$nu1)
  
  ## Longitudinal fixed effects
  long <- paste0("Y.", names(object$beta))
  
  ## Survival fixed effects
  survival1 <- paste0("T.", names(object$gamma1))
  
  ## Random-effect names used in association parameters
  re.names <- c("(Intercept)", random[seq_len(p1a - 1)])
  
  ## Association parameter names for event type 1
  nu1 <- paste0("T.asso:", re.names, "_1")
  
  ## Flexible covariance parameter names: Sig11, Sig22, Sig12, etc.
  sig <- c(
    paste0("Sig", seq_len(p1a), seq_len(p1a)),
    unlist(lapply(seq_len(p1a - 1), function(i) {
      paste0("Sig", i, (i + 1):p1a)
    }))
  )
  
  if (object$CompetingRisk) {
    survival2 <- paste0("T.", names(object$gamma2))
    nu2 <- paste0("T.asso:", re.names, "_2")
    
    par.names <- c(
      long,
      survival1,
      survival2,
      nu1,
      nu2,
      "Y.sigma^2",
      sig
    )
  } else {
    par.names <- c(
      long,
      survival1,
      nu1,
      "Y.sigma^2",
      sig
    )
  }
  
  colnames(vcov) <- par.names
  rownames(vcov) <- par.names
  
  vcov
}

