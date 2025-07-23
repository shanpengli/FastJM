##' @title Print mvjmcs
##' @name print
##' @aliases print.mvjmcs
##' @param x Object of class 'mvjmcs'.
##' @param digits the number of significant digits to use when printing.
##' @param ... Further arguments passed to or from other methods.
##' @return a summary of data, joint model, log likelihood, and parameter estimates.
##' @seealso \code{\link{mvjmcs}}
##' @export
print.mvjmcs <- function(x, digits = 4, ...) {
  if (!inherits(x, "mvjmcs"))
    stop("Use only with 'mvjmcs' objects.\n")
  
  cat("\nCall:\n", sprintf(format(paste(deparse(x$call, width.cutoff = 500), collapse = ""))), "\n\n")
  
  numBio <- length(x$sigma)
  
  if (x$CompetingRisk) {
    cat("Data Summary:\n")
    cat("Number of observations:", nrow(x$ydata), "\n")
    cat("Number of groups:", nrow(x$cdata), "\n\n")
    cat("Proportion of competing risks: \n")
    for (i in 1:2) {
      cat("Risk", i, ":", round(x$PropEventType[i+1, 2]/nrow(x$cdata)*100, 2), "%\n")
    }
    cat("\nModel Type: joint modeling of multivariate longitudinal continuous and competing risks data", "\n\n")
    cat("Model summary:\n")
    cat("Longitudinal process: linear mixed effects model\n")
    cat("Event process: cause-specific Cox proportional hazard model with non-parametric baseline hazard\n\n")
    
    cat("Fixed effects in the longitudinal sub-model: ",
        sprintf(format(paste(deparse(x$LongitudinalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    
    
    
    # ~~~~~~~~~~~~~~~~~~~~~~~
    # print beta coefficients
    # ~~~~~~~~~~~~~~~~~~~~~~~
    dat <- data.frame(x$beta, x$sebeta, x$beta/x$sebeta, 2 * pnorm(-abs(x$beta/x$sebeta)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    cat("\n")
    
    # ~~~~~~~~~~~~~~~~~
    # print sigma^2 est
    # ~~~~~~~~~~~~~~~~~
    dat <- data.frame(x$sigma, x$sesigma, x$sigma/x$sesigma, 2 * pnorm(-abs(x$sigma/x$sesigma)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    
    tempName <- c()
    for(g in 1:numBio){
      tempName[g] <- paste0("sigma^2_","bio", g)
    }
    
    rownames(dat) <- tempName
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    
    cat("\nFixed effects in the survival sub-model: ",
        sprintf(format(paste(deparse(x$SurvivalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    dat <- data.frame(x$gamma1, x$segamma1, x$gamma1/x$segamma1, 2 * pnorm(-abs(x$gamma1/x$segamma1)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    
    dat2 <- data.frame(x$gamma2, x$segamma2, x$gamma2/x$segamma2, 2 * pnorm(-abs(x$gamma2/x$segamma2)))
    colnames(dat2) <- c("Estimate", "SE", "Z value", "p-val")
    dat <- rbind(dat, dat2)
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    cat("\nAssociation parameters:                 \n")
    dat <- data.frame(x$alpha1, x$sealpha1, x$alpha1/x$sealpha1, 2 * pnorm(-abs(x$alpha1/x$sealpha1)))
    subdat <- data.frame(x$alpha2, x$sealpha2, x$alpha2/x$sealpha2, 2 * pnorm(-abs(x$alpha2/x$sealpha2)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    colnames(subdat) <- c("Estimate", "SE", "Z value", "p-val")
    dat <- rbind(dat, subdat)
    tempName <- c()
    ind = 1
    
    for(g in 1:numBio){
      pRE <- length(all.vars(x$random[[g]]))
      if (pRE == 1){
        tempName[ind] <- paste0("(Intercept)_1","bio",g)
      } else{
        temp <- c(paste0("(Intercept)_1","bio",g))
        tempName[ind:(ind+pRE-1)] <- c(temp, paste0(all.vars(x$random[[g]])[-pRE],"_1","bio", g))
      }  
      ind <- ind + pRE
    }
    
    for(g in 1:numBio){
      
      pRE <- length(all.vars(x$random[[g]]))
      if (pRE == 1){
        tempName[ind] <- paste0("(Intercept)_2","bio",g)
      } else{
        temp <- c(paste0("(Intercept)_2","bio",g))
        tempName[ind:(ind+pRE-1)] <- c(temp, paste0(all.vars(x$random[[g]])[-pRE],"_2","bio", g))
      }  
      ind <- ind + pRE
      
      
    }
    
    rownames(dat) <- tempName
    
    
    
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    cat("\n")
    
    cat("\nRandom effects:                 \n")
    for(g in 1:numBio){
      cat("  bio",g, ": ", format(as.formula(x$random[[g]])), "\n")
    }
    
    
    pREtotal <- length(diag(x$Sig))
    
    dat <- matrix(0, nrow = pREtotal * (pREtotal + 1)/2, ncol = 4)
    for(i in 1:pREtotal){
      dat[i,] <- c( x$Sig[i,i], x$seSig[i,i], x$Sig[i,i]/x$seSig[i,i], 2 * pnorm(-abs(x$Sig[i,i]/x$seSig[i,i])))
    }
    
    re <- c()
    pREvec <- c()
    for(g in 1:numBio){
      pREvec[g] <- length(all.vars(x$random[[g]]))
    }
    
    pREind <- 1
    
    
    for(g in 1:numBio){
      pRE <- pREvec[g]
      temp <- all.vars(x$random[[g]])
      ID <- all.vars(x$random[[g]])[length(all.vars(x$random[[g]]))]
      if (length(temp) == 1) {
        re[pREind:(pREind+pRE-1)] <- paste0("Intercept", g)
      } else {
        re[pREind:(pREind+pRE-1)] <- paste0(c("Intercept", 
                                              temp[-length(temp)]), g)
      }
      pREind <- pREind+pRE
    }
    
    name <- re
    
    ind <- pREtotal+1
    for (i in 1:(pREtotal-1)) {
      for (j in (i+1):(pREtotal)){
        dat[ind, ] <- c(x$Sig[i, j], x$seSig[i, j], x$Sig[i, j]/x$seSig[i, j],
                        2 * pnorm(-abs(x$Sig[i, j]/x$seSig[i, j])))
        name[ind] <- paste0(re[i], ":", re[j])
        ind <- ind + 1
      }
    }
    
    dat <- as.data.frame(dat)
    
    
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    rownames(dat) <- name
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
  } else {
    cat("Data Summary:\n")
    cat("Number of observations:", nrow(x$ydata), "\n")
    cat("Number of groups:", nrow(x$cdata), "\n\n")
    cat("Proportion of events:", round(x$PropEventType[2, 2]/nrow(x$cdata)*100, 2), "%\n")
    cat("\nModel Type: joint modeling of multivariate longitudinal continuous and survival data", "\n\n")
    cat("Model summary:\n")
    cat("Longitudinal process: linear mixed effects model\n")
    cat("Event process: Cox proportional hazard model with non-parametric baseline hazard\n\n")
    cat("Fixed effects in the longitudinal submodel: ",
        sprintf(format(paste(deparse(x$LongitudinalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    
    # ~~~~~~~~~~~~~~~~~~~~~~~
    # print beta coefficients
    # ~~~~~~~~~~~~~~~~~~~~~~~
    dat <- data.frame(x$beta, x$sebeta, x$beta/x$sebeta, 2 * pnorm(-abs(x$beta/x$sebeta)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    cat("\n")
    
    # ~~~~~~~~~~~~~~~~~
    # print sigma^2 est
    # ~~~~~~~~~~~~~~~~~
    dat <- data.frame(x$sigma, x$sesigma, x$sigma/x$sesigma, 2 * pnorm(-abs(x$sigma/x$sesigma)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    
    tempName <- c()
    for(g in 1:numBio){
      tempName[g] <- paste0("sigma^2_","bio", g)
    }
    
    rownames(dat) <- tempName
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    
    cat("\nFixed effects in the survival sub-model: ",
        sprintf(format(paste(deparse(x$SurvivalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    dat <- data.frame(x$gamma1, x$segamma1, x$gamma1/x$segamma1, 2 * pnorm(-abs(x$gamma1/x$segamma1)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    
    # ~~~~~
    # Alpha
    # ~~~~~
    cat("\nAssociation parameters:                 \n")
    dat <- data.frame(x$alpha1, x$sealpha1, x$alpha1/x$sealpha1, 2 * pnorm(-abs(x$alpha1/x$sealpha1)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    random <- all.vars(x$random)
    tempName <- c()
    ind = 1
    
    for(g in 1:numBio){
      pRE <- length(all.vars(x$random[[g]]))
      if (pRE == 1){
        tempName[ind] <- paste0("(Intercept)_1","bio",g)
      } else{
        temp <- c(paste0("(Intercept)_1","bio",g))
        tempName[ind:(ind+pRE-1)] <- c(temp, paste0(all.vars(x$random[[g]])[-pRE],"_1","bio", g))
      }  
      ind <- ind + pRE
    }
    
    rownames(dat) <- tempName
    
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    cat("\n")
    
    cat("\nRandom effects:                 \n")
    for(g in 1:numBio){
      cat("  bio",g, ": ", format(as.formula(x$random[[g]])), "\n")
    }
    
    
    pREtotal <- length(diag(x$Sig))
    
    dat <- matrix(0, nrow = pREtotal * (pREtotal + 1)/2, ncol = 4)
    for(i in 1:pREtotal){
      dat[i,] <- c( x$Sig[i,i], x$seSig[i,i], x$Sig[i,i]/x$seSig[i,i], 2 * pnorm(-abs(x$Sig[i,i]/x$seSig[i,i])))
    }
    
    re <- c()
    pREvec <- c()
    for(g in 1:numBio){
      pREvec[g] <- length(all.vars(x$random[[g]]))
    }
    
    pREind <- 1
    
    
    for(g in 1:numBio){
      pRE <- pREvec[g]
      temp <- all.vars(x$random[[g]])
      ID <- all.vars(x$random[[g]])[length(all.vars(x$random[[g]]))]
      if (length(temp) == 1) {
        re[pREind:(pREind+pRE-1)] <- paste0("Intercept", g)
      } else {
        re[pREind:(pREind+pRE-1)] <- paste0(c("Intercept", 
                                              temp[-length(temp)]), g)
      }
      pREind <- pREind+pRE
    }
    
    name <- re
    
    ind <- pREtotal+1
    for (i in 1:(pREtotal-1)) {
      for (j in (i+1):(pREtotal)){
        dat[ind, ] <- c(x$Sig[i, j], x$seSig[i, j], x$Sig[i, j]/x$seSig[i, j],
                        2 * pnorm(-abs(x$Sig[i, j]/x$seSig[i, j])))
        name[ind] <- paste0(re[i], ":", re[j])
        ind <- ind + 1
      }
    }
    
    dat <- as.data.frame(dat)
    
    
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    rownames(dat) <- name
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
  }
}
