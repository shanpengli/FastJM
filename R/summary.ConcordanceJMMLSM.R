##' @title Print ConcordanceJMMLSM
##' @name summary.ConcordanceJMMLSM
##' @aliases summary.ConcordanceJMMLSM
##' @param object object of class 'ConcordanceJMMLSM'.
##' @param digits number of digits of decimal to be printed. 
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{JMMLSM}, \link{survfitJMMLSM}}
##' @export
##' 

summary.ConcordanceJMMLSM <- function (object, digits = 4, ...) {
  if (!inherits(object, "ConcordanceJMMLSM"))
    stop("Use only with 'ConcordanceJMMLSM' xs.\n") 
  
  if (is.null(object$Concordance.cv)) {
    stop("The cross validation fails. Please try using a different seed number.")
  } else {
    
    if (length(object$Concordance.cv) == object$n.cv && sum(mapply(is.null, object$Concordance.cv)) == 0) {
      sum <- 0
      for (j in 1:object$n.cv) {
        sum <- sum + object$Concordance.cv[[j]]
      }
      sum <- sum/object$n.cv
      
      if (object$CompetingRisk) {
        Concordance1 <- sum[1]
        Concordance2 <- sum[2]
        ExpectedConcordance <- data.frame(Concordance1, Concordance2)
        colnames(ExpectedConcordance) <- c("Concordance1", "Concordance2")
      } else {
        Concordance <- sum[1]
        ExpectedConcordance <- data.frame(Concordance)
        colnames(ExpectedConcordance) <- c("Concordance")
      }
      return(ExpectedConcordance)
    } else {
      stop("The cross validation fails. Please try using a different seed number.")
    }
    
  }

}