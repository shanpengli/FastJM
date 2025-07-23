##' @title Anova Method for Fitted Joint Models
##' @name summary
##' @aliases summary.mvjmcs
##' @description Produce result summaries of a joint model fit. 
##' @param object an object inheriting from class \code{mvjmcs}.
##' @param process for which model (i.e., longitudinal model or survival model) to extract the estimated coefficients.
##' @param digits the number of significant digits to use when printing. Default is 4.
##' @param ... further arguments passed to or from other methods.
##' @return A table to summarize the model results.
##' @seealso \code{\link{mvjmcs}}
##' @export
##' 

summary.mvjmcs <- function(object, process = c("Longitudinal", "Event"), digits = 4, ...) {
  
  if (!inherits(object, "mvjmcs"))
    stop("Use only with 'mvjmcs' objects.\n")
  
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
    
    tempName <- c()
    numBio <- length(object$sigma)
    for(g in 1:numBio){
      tempName[g] <- paste0("sigma^2_","bio", g)
    }
    Estimate <- object$sigma
    SE <- object$sesigma
    LowerLimit <- Estimate - 1.96 * SE
    UpperLimit <- Estimate + 1.96 * SE
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    outsigma <- data.frame(Estimate, SE, LowerLimit, UpperLimit, pval)
    rownames(outsigma) <- tempName
    outsigma <- cbind(rownames(outsigma), outsigma)
    rownames(outsigma) <- NULL
    
    names(outsigma) <- c("Longitudinal", "coef", "SE", "95%Lower", "95%Upper", "p-values")
    
    outsigma[, 2:ncol(outsigma)] <- round(outsigma[, 2:ncol(outsigma)], digits = digits)
    outsigma[, ncol(outsigma)] <- format(outsigma[, ncol(outsigma)], scientific = FALSE)
    
    out <- rbind(out, outsigma)
    
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
    
    ##alpha
    Estimate <- object$alpha1
    if (length(Estimate) == 2) names(Estimate) <- c("alpha1_1", "alpha1_2")
    if (length(Estimate) == 1) names(Estimate) <- c("alpha1_1")
    SE <- object$sealpha1
    LowerLimit <- Estimate - 1.96 * SE
    expLL <- exp(LowerLimit)
    UpperLimit <- Estimate + 1.96 * SE
    expUL <- exp(UpperLimit)
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, 
                      expLL, expUL, pval)
    # out <- cbind(rownames(out), out)
    rownames(out) <- NULL
    colnames(out)[1] <- "Parameter"
    
    numBio <- length(object$sigma)
    tempName <- c()
    ind <- 1
    for(g in 1:numBio){
      pRE <- length(all.vars(object$random[[g]]))
      if (pRE == 1){
        tempName[ind] <- paste0("(Intercept)_1","bio",g)
      } else{
        temp <- c(paste0("(Intercept)_1","bio",g))
        tempName[ind:(ind+pRE-1)] <- c(temp, paste0(all.vars(object$random[[g]])[-pRE],"_1","bio", g))
      }  
      ind <- ind + pRE
    }
    
    if (object$CompetingRisk) {
      Estimate <- object$alpha2
      if(!is.null(Estimate)){
        for(g in 1:numBio){
          
          pRE <- length(all.vars(object$random[[g]]))
          if (pRE == 1){
            tempName[ind] <- paste0("(Intercept)_2","bio",g)
          } else{
            temp <- c(paste0("(Intercept)_2","bio",g))
            tempName[ind:(ind+pRE-1)] <- c(temp, paste0(all.vars(object$random[[g]])[-pRE],"_2","bio", g))
          }  
          ind <- ind + pRE
        }
        
        SE <- object$sealpha2
        LowerLimit <- Estimate - 1.96 * SE
        expLL <- exp(LowerLimit)
        UpperLimit <- Estimate + 1.96 * SE
        expUL <- exp(UpperLimit)
        zval = (Estimate/SE)
        pval = 2 * pnorm(-abs(zval))
        out2 <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, 
                           expLL, expUL, pval)
      }
      # print(out2)
      # out2 <- cbind(rownames(out2), out2)
      rownames(out2) <- NULL
      colnames(out2)[1] <- "Parameter"
      out <- rbind(out, out2)
    }
    out <- cbind(tempName, out)
    
    names(out) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", 
                    "95%exp(Lower)", "95%exp(Upper)", "p-values")
    
    out <- rbind(outgamma, out)
    
    out[, 2:ncol(out)] <- round(out[, 2:ncol(out)], digits = digits)
    out[, ncol(out)] <- format(out[, ncol(out)], scientific = FALSE)
    
    return(out)
  }
  
  
}
