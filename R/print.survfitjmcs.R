##' @title Print survfitjmcs
##' @name print.survfitjmcs
##' @aliases print.survfitjmcs
##' @param x x of class 'survfitjmcs'.
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}, \link{survfitjmcs}}
##' @export
##' 
print.survfitjmcs <- function (x, ...) {
  if (!inherits(x, "survfitjmcs"))
    stop("Use only with 'survfitjmcs' xs.\n")

  if (x$simulate) {
    cat("\nPrediction of Conditional Probabilities of Event\n\tbased on", x$M, "Monte Carlo samples\n\n")
    f <- function (d, t) {
      a <- matrix(1, nrow = 1, ncol = 5)
      a[1, 1] <- t 
      a <- as.data.frame(a)
      colnames(a) <- colnames(d)
      d <- rbind(a, d)
      d
    }
    
    f.CR <- function (d, t) {
      a <- matrix(0, nrow = 1, ncol = 5)
      a[1, 1] <- t 
      a <- as.data.frame(a)
      
      colnames(a) <- colnames(d[[1]])
      for (i in 1:2) {
        d[[i]] <- rbind(a, d[[i]])
      }
      d
    }
    x$Last.time <- as.data.frame(x$Last.time)
    if (!x$CompetingRisk) {
      print(mapply(f, x$Pred, x$Last.time[, 2], SIMPLIFY = FALSE))
    } else {
      print(mapply(f.CR, x$Pred, x$Last.time[, 2], SIMPLIFY = FALSE))
    }
  } else {
    if (!x$CompetingRisk) {
      cat("\nPrediction of Conditional Survival Probabilities\n", "(Confidence interval not available)\n")
      print(x$Pred)
    } else {
      cat("\nPrediction of Conditional Cumulative Incidence Rate\n", "(Confidence interval not available)\n")
      print(x$Pred)
    }

  }
  invisible(x)
  }