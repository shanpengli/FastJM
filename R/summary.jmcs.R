##' @title Summarize fitted joint models
##' @name summary
##' @description Summarize parameter estimates for the longitudinal
##'   and event submodels of a fitted joint model.
##' @param object A fitted model of class \code{jmcs},
##'   \code{JMMLSM}, or \code{mvjmcs}.
##' @param digits Number of decimal places used for rounding
##'   numeric results. Default is \code{4}.
##' @param ... Further arguments passed to or from other methods.
##' @return A named list containing \code{Longitudinal} and
##'   \code{Event} data frames with parameter estimates and
##'   associated uncertainty, returned invisibly after printing.
##' @seealso \code{\link{jmcs}}, \code{\link{JMMLSM}},
##'   \code{\link{mvjmcs}}
##' @export

summary.jmcs <- function(object, digits = 4, ...) {
  
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  
  ##Estimates of betas
  Estimate <- object$beta
  SE <- object$sebeta
  LowerLimit <- Estimate - 1.96 * SE
  UpperLimit <- Estimate + 1.96 * SE
  zval = (Estimate/SE)
  pval = 2 * pnorm(-abs(zval))
  out <- data.frame(Estimate, SE, LowerLimit, UpperLimit, pval)
  out <- cbind(rownames(out), out)
  rownames(out) <- NULL
  
  names(out) <- c("Longitudinal", "coef", "SE", "95%Lower", "95%Upper", "p-values")
  
  out[, 2:ncol(out)] <- round(out[, 2:ncol(out)], digits = digits)
  out[, ncol(out)] <- format(out[, ncol(out)], scientific = FALSE)
  
  Longitudinal <- out
  cat("\nLongitudinal submodel:\n")
  print(Longitudinal, row.names = FALSE)
  
  ##gamma
  Estimate <- object$gamma1
  SE <- object$segamma1
  LowerLimit <- Estimate - 1.96 * SE
  expLL <- exp(LowerLimit)
  UpperLimit <- Estimate + 1.96 * SE
  expUL <- exp(UpperLimit)
  zval = (Estimate/SE)
  pval = 2 * pnorm(-abs(zval))
  out <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, 
                    expLL, expUL, pval)
  out <- cbind(rownames(out), out)
  rownames(out) <- NULL
  colnames(out)[1] <- "Parameter"
  outgamma <- out
  
  if (object$CompetingRisk) {
    Estimate <- object$gamma2
    SE <- object$segamma2
    LowerLimit <- Estimate - 1.96 * SE
    expLL <- exp(LowerLimit)
    UpperLimit <- Estimate + 1.96 * SE
    expUL <- exp(UpperLimit)
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out2 <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, 
                       expLL, expUL, pval)
    out2 <- cbind(rownames(out2), out2)
    rownames(out2) <- NULL
    colnames(out2)[1] <- "Parameter"
    outgamma <- rbind(out, out2)
  }
  names(outgamma) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", 
                       "95%exp(Lower)", "95%exp(Upper)", "p-values")
  
  ##nu
  Estimate <- as.numeric(object$nu1)
  SE       <- as.numeric(object$senu1)
  
  stopifnot(length(Estimate) == length(SE))
  
  # Ensure names: keep existing if present; otherwise assign nu1_1, nu1_2, ...
  if (is.null(names(Estimate)) || any(names(Estimate) == "")) {
    names(Estimate) <- paste0("nu1_", seq_along(Estimate))
  }
  
  LowerLimit <- Estimate - 1.96 * SE
  expLL <- exp(LowerLimit)
  UpperLimit <- Estimate + 1.96 * SE
  expUL <- exp(UpperLimit)
  zval = (Estimate/SE)
  pval = 2 * pnorm(-abs(zval))
  out <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, 
                    expLL, expUL, pval)
  out <- cbind(rownames(out), out)
  rownames(out) <- NULL
  colnames(out)[1] <- "Parameter"
  outnu <- out
  
  if(object$CompetingRisk) {
    
    Estimate <- as.numeric(object$nu2)
    SE       <- as.numeric(object$senu2)
    
    stopifnot(length(Estimate) == length(SE))
    
    # Ensure names: keep existing if present; otherwise assign nu1_1, nu1_2, ...
    if (is.null(names(Estimate)) || any(names(Estimate) == "")) {
      names(Estimate) <- paste0("nu2_", seq_along(Estimate))
    }
    
    LowerLimit <- Estimate - 1.96 * SE
    expLL <- exp(LowerLimit)
    UpperLimit <- Estimate + 1.96 * SE
    expUL <- exp(UpperLimit)
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out2 <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, 
                       expLL, expUL, pval)
    out2 <- cbind(rownames(out2), out2)
    rownames(out2) <- NULL
    colnames(out2)[1] <- "Parameter"
    outnu <- rbind(out, out2)
  }
  names(outnu) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", 
                    "95%exp(Lower)", "95%exp(Upper)", "p-values")
  
  out <- rbind(outgamma, outnu)
  
  out[, 2:ncol(out)] <- round(out[, 2:ncol(out)], digits = digits)
  out[, ncol(out)] <- format(out[, ncol(out)], scientific = FALSE)
  
  Event <- out
  cat("\nEvent submodel:\n")
  print(Event, row.names = FALSE)
  
  invisible(list(
    Longitudinal = Longitudinal,
    Event = Event
  ))
  
}