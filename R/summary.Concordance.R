##' @title Summarize Concordance
##' @name summary.Concordance
##' @aliases summary.Concordance
##' @method summary Concordance
##'
##' @param object An object of class 'Concordance'.
##' @param digits Number of decimal places to be printed. Default is \code{4}.
##' @param ... Further arguments passed to or from other methods.
##'
##' @return A data frame containing the cross-validation averaged concordance estimate.
##'
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{Concordance}}, \code{\link{jmcs}}, \code{\link{JMMLSM}}
##'
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