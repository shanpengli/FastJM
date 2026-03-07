##' @title Print PEJMMLSM
##' @name summary.PEJMMLSM
##' @aliases summary.PEJMMLSM
##' @param object object of class 'PEJMMLSM'.
##' @param error a character string that specifies the loss function. 
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{JMMLSM}, \link{survfitJMMLSM}}
##' @export
##' 

summary.PEJMMLSM <- function (object, error = c("MAE", "Brier"), ...) {
  if (!inherits(object, "PEJMMLSM"))
    stop("Use only with 'PEJMMLSM' xs.\n") 
  if (!error %in% c("MAE", "Brier"))
    stop("Please choose one of the following options: MAE or Brier.")
  
  if (error == "Brier") {
    if (is.null(object$Brier.cv)) {
      stop("The cross validation fails. Please try using a different seed number.")
    } else {
      
      if (length(object$Brier.cv) == object$n.cv && sum(mapply(is.null, object$Brier.cv)) == 0) {
        sum <- 0
        for (j in 1:object$n.cv) {
          sum <- sum + object$Brier.cv[[j]]
        }
        sum <- sum/object$n.cv
        
        if (object$CompetingRisk) {
          mean.Brier1 <- sum[, 1]
          mean.Brier2 <- sum[, 2]
          ExpectedBrier <- data.frame(object$horizon.time, mean.Brier1, mean.Brier2)
          colnames(ExpectedBrier) <- c("Horizon Time", "Brier Score 1", "Brier Score 2")
        } else {
          mean.Brier1 <- sum[, 1]
          ExpectedBrier <- data.frame(object$horizon.time, mean.Brier1)
          colnames(ExpectedBrier) <- c("Horizon Time", "Brier Score")
        }
        cat("\nExpected Brier Score at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
        print(ExpectedBrier)
      } else {
        stop("The cross validation fails. Please try using a different seed number.")
      }
      
    }
  } else {
    if (is.null(object$MAE.cv)) {
      stop("The cross validation fails. Please try using a different seed number.")
    } else {
      
      if (length(object$MAE.cv) == object$n.cv && sum(mapply(is.null, object$MAE.cv)) == 0) {
        sum <- 0
        for (j in 1:object$n.cv) {
          sum <- sum + object$MAE.cv[[j]]
        }
        sum <- sum/object$n.cv
        if (object$CompetingRisk) {
          mean.MAE1 <- sum[, 1]
          mean.MAE2 <- sum[, 2]
          
          ExpectedAE <- data.frame(object$horizon.time, mean.MAE1, mean.MAE2)
          colnames(ExpectedAE) <- c("Horizon Time", "MAE 1", "MAE 2")
        } else {
          mean.MAE1 <- sum[, 1]
          
          ExpectedAE <- data.frame(object$horizon.time, mean.MAE1)
          colnames(ExpectedAE) <- c("Horizon Time", "MAE")
        }
        cat("\nExpected mean absolute error at the landmark time of", object$landmark.time, "\nbased on", object$n.cv, "fold cross validation\n")
        print(ExpectedAE)
      } else {
        stop("The cross validation fails. Please try using a different seed number.")
      }
      
    }
    
  }
  
  
}