##' @title Print AUCjmcs
##' @name summary.AUCjmcs
##' @aliases summary.AUCjmcs
##' @param object object of class 'AUCjmcs'.
##' @param digits number of digits of decimal to be printed. 
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}, \link{survfitjmcs}}
##' @export
##' 

summary.AUCjmcs <- function (object, digits = 4, ...) {
  if (!inherits(object, "AUCjmcs"))
    stop("Use only with 'AUCjmcs' xs.\n") 
  
  if (is.null(object$AUC.cv)) {
    stop("The cross validation fails. Please try using a different seed number.")
  } else {
    
    if (length(object$AUC.cv) == object$n.cv && sum(mapply(is.null, object$AUC.cv)) == 0) {
      sum <- 0
      for (j in 1:object$n.cv) {
        sum <- sum + object$AUC.cv[[j]]
      }
      sum <- sum/object$n.cv
      
      if (object$CompetingRisk) {
        AUC1 <- sum[, 1]
        AUC2 <- sum[, 2]
        ExpectedAUC <- data.frame(object$horizon.time, AUC1, AUC2)
        colnames(ExpectedAUC) <- c("Horizon Time", "AUC1", "AUC2")
      } else {
        AUC <- sum[, 1]
        ExpectedAUC <- data.frame(object$horizon.time, AUC)
        colnames(ExpectedAUC) <- c("Horizon Time", "AUC")
      }
      cat("\nExpected AUC at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
      return(ExpectedAUC)
    } else {
      stop("The cross validation fails. Please try using a different seed number.")
    }
    
  }
  
  
  
  
}