##' @title Print MAEQjmcs
##' @name summary.MAEQjmcs
##' @aliases summary.MAEQjmcs
##' @param x x of class 'MAEQjmcs'.
##' @param digits number of decimal points to be rounded.
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}, \link{survfitjmcs}}
##' @export
##' 

summary.MAEQjmcs <- function (x, digits = 3, ...) {
  if (!inherits(x, "MAEQjmcs"))
    stop("Use only with 'MAEQjmcs' xs.\n") 
  
  if (is.null(x$MAEQ.cv)) {
    stop("The cross validation fails. Please try using a different seed number.")
  } else {
    if(length(x$MAEQ.cv) == x$n.cv && sum(mapply(is.null, x$MAEQ.cv)) == 0) {
      if (x$CompetingRisk) {
        sum <- as.data.frame(matrix(0, nrow = length(x$horizon.time), ncol = 3))
        sum[, 1] <- x$horizon.time
        colnames(sum) <- c("Horizon Time", "CIF1", "CIF2")
        for (i in 1:length(x$horizon.time)) {
          for (j in 1:x$n.cv) {
            sum[i, 2] <- sum[i, 2] + sum(abs(x$MAEQ.cv[[j]]$AllCIF1[[i]][, 1] - 
                                               x$MAEQ.cv[[j]]$AllCIF1[[i]][, 2])) 
            sum[i, 3] <- sum[i, 3] + sum(abs(x$MAEQ.cv[[j]]$AllCIF2[[i]][, 1] - 
                                               x$MAEQ.cv[[j]]$AllCIF2[[i]][, 2])) 
          }
        }
        sum[, -1] <- sum[, -1]/x$n.cv
      } else {
        sum <- as.data.frame(matrix(0, nrow = length(x$horizon.time), ncol = 2))
        sum[, 1] <- x$horizon.time
        colnames(sum) <- c("Horizon Time", "SurvProb")
        for (i in 1:length(x$horizon.time)) {
          for (j in 1:x$n.cv) {
            sum[i, 2] <- sum[i, 2] + sum(abs(x$MAEQ.cv[[j]]$AllSurv[[i]][, 1] - 
                                               x$MAEQ.cv[[j]]$AllSurv[[i]][, 2])) 
          }
        }
        sum[, -1] <- sum[, -1]/x$n.cv
      }
      sum[, -1] <- round(sum[, -1], digits)
      cat("\nSum of absolute error across quintiles of predicted risk scores at the landmark time of", x$landmark.time, "\nbased on", x$n.cv, "fold cross validation\n")
      print(sum)
      invisible(x)
    } else {
      stop("The cross validation fails. Please try using a different seed number.")
    }
  }
  
}