##' @title Print PAjmcs
##' @name summary.PAjmcs
##' @aliases summary.PAjmcs
##' @param object object of class 'PAjmcs'.
##' @param digits number of decimal points to be rounded.
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}, \link{survfitjmcs}}
##' @export
##' 

summary.PAjmcs <- function (object, digits = 3, ...) {
  if (!inherits(object, "PAjmcs"))
    stop("Use only with 'PAjmcs' xs.\n") 

  
  if (is.null(object$PA.cv)) {
    stop("The cross validation fails. Please try using a different seed number.")
  } else {
    
    if (length(object$PA.cv) == object$n.cv && sum(mapply(is.null, object$PA.cv)) == 0) {
      sum <- matrix(0, nrow = length(object$horizon.time), ncol = 2)
      for (j in 1:object$n.cv) {
        sum <- sum + object$PA.cv[[j]]
      }
      sum <- sum/object$n.cv
      
      if (object$CompetingRisk) {
        mean.PA1 <- sum[, 1]
        mean.PA2 <- sum[, 2]
        sum <- round(sum, digits)
        ExpectedPA <- data.frame(object$horizon.time, mean.PA1, mean.PA2)
        colnames(ExpectedPA) <- c("Horizon Time", "psuedo R square 1", "psuedo R square 2")
      } else {
        mean.PA1 <- sum[, 1]
        sum <- round(sum, digits)
        ExpectedPA <- data.frame(object$horizon.time, mean.PA1)
        colnames(ExpectedPA) <- c("Horizon Time", "psuedo R square")
      }
      cat("\nExpected psuedo R square at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
      print(ExpectedPA)
    } else {
      stop("The cross validation fails. Please try using a different seed number.")
    }
    
  }
  
  
  
}