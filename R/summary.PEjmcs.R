##' @title Print PEjmcs
##' @name summary.PEjmcs
##' @aliases summary.PEjmcs
##' @param x x of class 'PEjmcs'.
##' @param error a character string that specifies the loss function. 
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}, \link{survfitjmcs}}
##' @export
##' 

summary.PEjmcs <- function (x, error = c("MAE", "Brier"), ...) {
  if (!inherits(x, "PEjmcs"))
    stop("Use only with 'PEjmcs' xs.\n") 
  if (!error %in% c("MAE", "Brier"))
    stop("Please choose one of the following options: MAE or Brier.")
  
  if (error == "Brier") {
    if (is.null(x$Brier.cv)) {
      stop("The cross validation fails. Please try using a different seed number.")
    } else {
      
      if (length(x$Brier.cv) == x$n.cv && sum(mapply(is.null, x$Brier.cv)) == 0) {
        sum <- 0
        for (j in 1:x$n.cv) {
          sum <- sum + x$Brier.cv[[j]]
        }
        sum <- sum/x$n.cv
        mean.Brier1 <- sum[, 1]
        mean.Brier2 <- sum[, 2]
        
        ExpectedBrier <- data.frame(x$horizon.time, mean.Brier1, mean.Brier2)
        colnames(ExpectedBrier) <- c("Horizon Time", "Failure 1", "Failure 2")
        cat("\nExpected Brier Score at the landmark time of", x$landmark.time, "\nbased on", x$n.cv, "fold cross validation\n")
        print(ExpectedBrier)
        invisible(x)
      } else {
        stop("The cross validation fails. Please try using a different seed number.")
      }
      
    }
  } else {
    if (is.null(x$MAE.cv)) {
      stop("The cross validation fails. Please try using a different seed number.")
    } else {
      
      if (length(x$MAE.cv) == x$n.cv && sum(mapply(is.null, x$MAE.cv)) == 0) {
        sum <- 0
        for (j in 1:x$n.cv) {
          sum <- sum + x$MAE.cv[[j]]
        }
        sum <- sum/x$n.cv
        mean.MAE1 <- sum[, 1]
        mean.MAE2 <- sum[, 2]
        
        ExpectedAE <- data.frame(x$horizon.time, mean.MAE1, mean.MAE2)
        colnames(ExpectedAE) <- c("Horizon Time", "Failure 1", "Failure 2")
        cat("\nExpected mean absolute error at the landmark time of", x$landmark.time, "\nbased on", x$n.cv, "fold cross validation\n")
        print(ExpectedAE)
        invisible(x)
      } else {
        stop("The cross validation fails. Please try using a different seed number.")
      }
      
    }
    
  }
  
  
}