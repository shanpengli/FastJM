##' @title Fitted values for joint models
##' @name plot.jmcs
##' @aliases plot.jmcs
##' @description Plot Diagnostics for Joint Models.
##' @param x x of class 'jmcs'.
##' @param add.smooth logical; if \code{TRUE} a smooth line is superimposed in the "Residuals vs Fitted" plot.
##' @param ... further arguments passed to or from other methods.
##' @return The first two plots are longitudinal sub-model diagnostics and the last two are marginal survival function and marginal cumulative hazard.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @examples 
##' \donttest{
##' fit <- jmcs(ydata = ydata, cdata = cdata, 
##'             long.formula = response ~ time + gender + x1 + race, 
##'             surv.formula = Surv(surv, failure_type) ~ x1 + gender + x2 + race, 
##'             random =  ~ time| ID)
##' 
##' plot(fit)
##' }
##' 
##' @export
##' 

plot.jmcs <- function(x, add.smooth = getOption("add.smooth"), ...) {
  
  if (!inherits(x, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  
  op <- par(mfrow = c(2, 2))
  on.exit(par(op), add = TRUE)
  
  residuals <- x$fitted$resid
  fitted <- x$fitted$fitted
  
  plot(
    fitted, residuals,
    xlab = "Fitted Values",
    ylab = "Residuals",
    main = "Residuals vs Fitted"
  )
  if (isTRUE(add.smooth)) {
    abline(h = 0, lty = 3, col = "grey", lwd = 2)
    panel.smooth(fitted, residuals, lwd = 2)
  }
  
  qqnorm(
    residuals,
    ylab = "Standardized Residuals",
    main = "Normal Q-Q",
    ...
  )
  qqline(residuals, lty = 3, col = "grey50")
  
  marsurv <- as.data.frame(x$fittedSurv)
  
  plot(
    marsurv$V2 ~ marsurv$V1,
    type = "l",
    main = "Marginal Survival",
    ylab = "Survival Probability",
    xlab = "Time"
  )
  
  plot(
    -log(marsurv$V2) ~ marsurv$V1,
    type = "l",
    main = "Marginal Cumulative Hazard",
    ylab = "Cumulative Hazard",
    xlab = "Time"
  )
  
  invisible()
}

