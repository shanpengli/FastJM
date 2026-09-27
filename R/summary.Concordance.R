##' @title Summarize concordance estimates
##' @name summary.Concordance
##' @description Summarize concordance estimates averaged across
##'   cross-validation folds.
##' @param object An object of class \code{Concordance}.
##' @param digits Number of decimal places used for rounding.
##'   Default is \code{4}.
##' @param ... Further arguments passed to or from other methods.
##' @return A data frame containing the averaged concordance estimates.
##' @seealso \code{\link{Concordance}}
##' @export

summary.Concordance <- function(object, digits = 4, ...) {
  
  if (!inherits(object, "Concordance")) {
    stop("Use only with 'Concordance' objects.\n")
  }
  
  if (is.null(object$Concordance.cv)) {
    stop("The cross-validation fails. Please try using a different seed number.")
  }
  
  if (length(object$Concordance.cv) != object$n.cv ||
      any(vapply(object$Concordance.cv, is.null, logical(1)))) {
    stop("The cross-validation fails. Please try using a different seed number.")
  }
  
  Concordance.mean <- Reduce("+", object$Concordance.cv) / object$n.cv
  Concordance.mean <- round(Concordance.mean, digits = digits)
  
  if (object$CompetingRisk) {
    ExpectedConcordance <- data.frame(
      Concordance1 = Concordance.mean[1],
      Concordance2 = Concordance.mean[2]
    )
  } else {
    ExpectedConcordance <- data.frame(
      Concordance = Concordance.mean[1]
    )
  }
  
  return(ExpectedConcordance)
}