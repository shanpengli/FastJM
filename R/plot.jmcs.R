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
  
  if (!inherits(x, "jmcs")) {
    stop("Use only with 'jmcs' objects.\n")
  }
  
  op <- graphics::par(mfrow = c(2, 2))
  on.exit(graphics::par(op), add = TRUE)
  
  resid.values <- x$fitted$resid
  fitted <- x$fitted$fitted
  
  graphics::plot(
    fitted,
    resid.values,
    xlab = "Fitted Values",
    ylab = "Residuals",
    main = "Residuals vs Fitted"
  )
  
  if (isTRUE(add.smooth)) {
    graphics::abline(
      h = 0,
      lty = 3,
      col = "grey",
      lwd = 2
    )
    
    graphics::panel.smooth(
      fitted,
      resid.values,
      lwd = 2
    )
  }
  
  stats::qqnorm(
    resid.values,
    ylab = "Standardized Residuals",
    main = "Normal Q-Q",
    ...
  )
  
  stats::qqline(
    resid.values,
    lty = 3,
    col = "grey50"
  )
  
  marsurv <- as.data.frame(x$fittedSurv)
  
  cdata <- x$cdata
  surv.formula <- x$SurvivalSubmodel
  surv.var <- all.vars(surv.formula)
  
  survdata <- cdata[, surv.var[1:2], drop = FALSE]
  colnames(survdata) <- c("time", "status")
  
  km_fit <- survival::survfit(
    survival::Surv(time, status != 0) ~ 1,
    data = survdata
  )
  
  graphics::plot(
    marsurv$V2 ~ marsurv$V1,
    type = "l",
    lty = 1,
    lwd = 2,
    main = "Marginal Survival",
    ylab = "Survival Probability",
    xlab = "Time",
    ylim = c(0, 1)
  )
  
  graphics::lines(
    km_fit,
    lty = 2,
    lwd = 2,
    conf.int = FALSE
  )
  
  graphics::legend(
    "topright",
    legend = c("Estimated", "Kaplan-Meier"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
  
  na_fit <- survival::survfit(
    survival::Surv(time, status != 0) ~ 1,
    data = survdata,
    type = "fh2"
  )
  
  graphics::plot(
    marsurv$V1,
    -log(marsurv$V2),
    type = "l",
    lwd = 2,
    lty = 1,
    ylim = c(
      0,
      max(
        -log(marsurv$V2),
        na_fit$cumhaz,
        na.rm = TRUE
      )
    ),
    main = "Marginal Cumulative Hazard",
    ylab = "Cumulative Hazard",
    xlab = "Time"
  )
  
  graphics::lines(
    na_fit$time,
    na_fit$cumhaz,
    lty = 2,
    lwd = 2
  )
  
  graphics::legend(
    "topleft",
    legend = c("Estimated", "Nelson-Aalen"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
  
  invisible(x)
}