##' @title Random effects estimates for joint models
##' @name ranef
##' @description Extracts the posterior mean of the random effects for a fitted joint model.
##' @param object an object inheriting from class \code{jmcs} or \code{mvjmcs}.
##' @param ... further arguments passed to or from other methods.
##' @return a matrix of random effects estimates.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}}
##' @examples 
##' \donttest{
##' # a joint model fit
##' fit <- jmcs(ydata = ydata, cdata = cdata, 
##'             long.formula = response ~ time + gender + x1 + race, 
##'             surv.formula = Surv(surv, failure_type) ~ x1 + gender + x2 + race, 
##'             random =  ~ time| ID)
##' 
##' # extract random effects estimates
##' head(ranef(fit))
##' }
##' @export
##' 

ranef <- function(object, ...) {
  if (!inherits(object, "jmcs") && !inherits(object, "mvjmcs"))
    stop("Use only with 'mvjmcs' or jmcs' objects.\n")
  
  if(inherits(object, "jmcs")){
    rand <- all.vars(object$random)
    dat <- data.frame(t(object$FUNB))
    if (length(rand) == 1) {
      colnames(dat) <- c("(Intercept)") 
    } else {
      colnames(dat) <- c("(Intercept)", rand[-length(rand)]) 
    }
  }else{
    rand <- c()
    ind <- 1
    tempName <- c()
    numBio <- length(object$random)
    for(g in 1:numBio){
      pRE <- length(all.vars(object$random[[g]]))
      if (pRE == 1){
        tempName[ind] <- paste0("(Intercept)_","bio",g)
      } else{
        temp <- c(paste0("(Intercept)_","bio",g))
        tempName[ind:(ind+pRE-1)] <- c(temp, paste0(all.vars(object$random[[g]])[-pRE],"_","bio", g))
      }  
      ind <- ind + pRE
    }
    
    dat <- data.frame(matrix(unlist(object$pos.mode), 
                             nrow = length(object$pos.mode), byrow = TRUE))
    colnames(dat) <- tempName
  
  }

  return(dat)
}
