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
  
  f <- function (d, t) {
    a <- matrix(1, nrow = 1, ncol = 2)
    a[1, 1] <- t 
    a <- as.data.frame(a)
    colnames(a) <- colnames(d)
    d <- rbind(a, d)
    d
  }
  
  f.CR <- function (d, t) {
    a <- matrix(0, nrow = 1, ncol = 3)
    a[1, 1] <- t 
    a <- as.data.frame(a)
    colnames(a) <- colnames(d)
    d <- rbind(a, d)
    d
  }
  if (!is.null(x$quadpoint)) {
    cat("\nPrediction of Conditional Probabilities of Event\nbased on the pseudo-adaptive Guass-Hermite quadrature rule with", x$quadpoint,
        "quadrature points\n")
  } else {
    cat("\nPrediction of Conditional Probabilities of Event\nbased on the first order approximation\n")
  }
  if (!x$CompetingRisk) {
    print(mapply(f, x$Pred, x$Last.time[, 2], SIMPLIFY = FALSE))
  } else {
    print(mapply(f.CR, x$Pred, x$Last.time[, 2], SIMPLIFY = FALSE))
  }
  invisible(x)
  
}