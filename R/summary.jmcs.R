##' @title Anova Method for Fitted Joint Models
##' @name summary
##' @aliases summary.jmcs
##' @description Produce result summaries of a joint model fit. 
##' @param object an object inheriting from class \code{jmcs}.
##' @param process for which model (i.e., longitudinal model or survival model) to extract the estimated coefficients.
##' @param digits the number of significant digits to use when printing. Default is 4.
##' @param ... further arguments passed to or from other methods.
##' @return A table to summarize the model results.
##' @seealso \code{\link{jmcs}}
##' @export
##' 

summary.jmcs <- function(object, process = c("Longitudinal", "Event"), digits = 4, ...) {
  
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  
  if (process == "Longitudinal") {
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
    
    return(out)
    
  } else if (process == "Event") {
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
    names(outgamma) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", 
                         "95%exp(Lower)", "95%exp(Upper)", "p-values")
    
    ##nu
    Estimate <- object$nu1
    if (length(Estimate) == 2) names(Estimate) <- c("nu1_1", "nu1_2")
    if (length(Estimate) == 1) names(Estimate) <- c("nu1_1")
    SE <- object$senu1
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
    
    Estimate <- object$nu2
    if (length(Estimate) == 2) names(Estimate) <- c("nu2_1", "nu2_2")
    if (length(Estimate) == 1) names(Estimate) <- c("nu2_1")
    SE <- object$senu2
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
    names(outnu) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", 
                      "95%exp(Lower)", "95%exp(Upper)", "p-values")
    
    out <- rbind(outgamma, outnu)
    
    out[, 2:ncol(out)] <- round(out[, 2:ncol(out)], digits = digits)
    out[, ncol(out)] <- format(out[, ncol(out)], scientific = FALSE)
    
    return(out)
  }
  
  
}