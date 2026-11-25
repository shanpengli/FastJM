##' @title Print dynamic predictions from joint models
##' @name print.survfit
##' @description Print subject-specific conditional survival
##'   probabilities for single-failure models or conditional
##'   cumulative incidence probabilities for competing-risks models.
##' @param x A prediction object of class \code{survfitjmcs},
##'   \code{survfitJMMLSM}, or \code{survfitmvjmcs}.
##' @param ... Further arguments passed to or from other methods.
##' @return The input object \code{x}, invisibly.
##' @seealso \code{\link{survfitJM}}
##' @export

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
    cat("\nPrediction of Conditional Probabilities of Event\nbased on the pseudo-adaptive Gauss-Hermite quadrature rule with", x$quadpoint,
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