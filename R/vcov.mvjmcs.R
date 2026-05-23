##' @title Variance-covariance matrix of the estimated parameters for joint models
##' @name vcov
##' @aliases vcov.mvjmcs
##' @description Extract variance-covariance matrix for joint models.
##' @param object an object inheriting from class \code{mvjmcs}.
##' @param ... further arguments passed to or from other methods.
##' @return a matrix of variance covariance of all parameter estimates.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{mvjmcs}}
##' @export
##' 

vcov.mvjmcs <- function(object, ...) {
  if (!inherits(object, "mvjmcs"))
    stop("Use only with 'mvjmcs' objects.\n")

  vcov <- as.data.frame(object$vcov)
  numBio <- length(object$LongitudinalSubmodel)
  
  ## Longitudinal fixed effects
  long <- paste0("Y.", names(object$beta))
  
  sigma <- c()
  for(g in 1:numBio){
    sigma[g] <- paste0("Y.sigma^2_","bio", g)
  }
  
  ## Survival fixed effects
  survival1 <- paste0("T.", names(object$gamma1))
  
  ## Random-effect names used in association parameters
  re.names <- c()
  ind = 1
  for(g in 1:numBio){
    pRE <- length(all.vars(object$random[[g]]))
    if (pRE == 1){
      re.names[ind] <- paste0("(Intercept)_","bio",g)
    } else{
      temp <- c(paste0("(Intercept)_","bio",g))
      re.names[ind:(ind+pRE-1)] <- c(temp, paste0(all.vars(object$random[[g]])[-pRE],"_","bio", g))
    }  
    ind <- ind + pRE
  }
  
  ## Association parameter names for event type 1
  alpha1 <- paste0("T.asso:", re.names, "_1")
  
  ## Flexible covariance parameter names
  pREtotal <- nrow(object$Sig)
  
  re <- character(0)
  
  for (g in seq_len(numBio)) {
    temp <- all.vars(object$random[[g]])
    random.vars <- temp[-length(temp)]  # remove ID variable
    
    if (length(random.vars) == 0) {
      re <- c(re, paste0("Intercept", g))
    } else {
      re <- c(re, paste0(c("Intercept", random.vars), g))
    }
  }
  
  ## diagonal terms first
  sig <- paste0(re, ":", re)
  
  ## off-diagonal terms by lag order
  for (lag in seq_len(pREtotal - 1)) {
    for (i in seq_len(pREtotal - lag)) {
      j <- i + lag
      sig <- c(sig, paste0(re[i], ":", re[j]))
    }
  }
  
  if (object$CompetingRisk) {
    survival2 <- paste0("T.", names(object$gamma2))
    alpha2 <- paste0("T.asso:", re.names, "_2")
    
    par.names <- c(
      long,
      sigma,
      survival1,
      survival2,
      alpha1,
      alpha2,
      sig
    )
  } else {
    par.names <- c(
      long,
      sigma,
      survival1,
      alpha1,
      sig
    )
  }
  
  colnames(vcov) <- par.names
  rownames(vcov) <- par.names
  
  vcov
}

