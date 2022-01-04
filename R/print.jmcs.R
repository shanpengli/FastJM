##' @title Print jmcs
##' @name print
##' @aliases print.jmcs
##' @param x Object of class 'jmcs'.
##' @param digits the number of significant digits to use when printing. 
##' @param ... Further arguments passed to or from other methods.
##' @return a summary of data, joint model, log likelihood, and parameter estimates.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}}
##' @export
##' 
print.jmcs <- function(x, digits = 4, ...) {
  if (!inherits(x, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  
  cat("\nCall:\n", sprintf(format(paste(deparse(x$call, width.cutoff = 500), collapse = ""))), "\n\n")
  
  if (x$CompetingRisk) {
    cat("Data Summary:\n")
    cat("Number of observations:", nrow(x$ydata), "\n")
    cat("Number of groups:", nrow(x$cdata), "\n\n")
    cat("Proportion of competing risks: \n")
    for (i in 1:2) {
      cat("Risk", i, ":", round(x$PropEventType[i+1, 2]/nrow(x$cdata)*100, 2), "%\n")
    }
    cat("\nNumerical intergration:\n")
    cat(paste("Method:", x$Quad.method, "Guass-Hermite quadrature\n"))
    cat("Number of quadrature points: ", x$quadpoint, "\n")
    cat("\nModel Type: joint modeling of longitudinal continuous and competing risks data", "\n\n")
    cat("Model summary:\n")
    cat("Longitudinal process: linear mixed effects model\n")
    cat("Event process: cause-specific Cox proportional hazard model with non-parametric baseline hazard\n\n")
    cat("Loglikelihood: ", x$loglike, "\n\n")
    cat("Fixed effects in the longitudinal sub-model: ",
        sprintf(format(paste(deparse(x$LongitudinalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    
    dat <- data.frame(x$beta, x$sebeta, x$beta/x$sebeta, 2 * pnorm(-abs(x$beta/x$sebeta)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    cat("\n")
    dat <- data.frame(x$sigma, x$sesigma, x$sigma/x$sesigma, 2 * pnorm(-abs(x$sigma/x$sesigma)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    rownames(dat) <- "sigma^2"
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
    dat <- data.frame(x$nu1, x$senu1, x$nu1/x$senu1, 2 * pnorm(-abs(x$nu1/x$senu1)))
    subdat <- data.frame(x$nu2, x$senu2, x$nu2/x$senu2, 2 * pnorm(-abs(x$nu2/x$senu2)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    colnames(subdat) <- c("Estimate", "SE", "Z value", "p-val")
    dat <- rbind(dat, subdat)
    if (length(x$nu1) == 1) rownames(dat) <- c("nu1_1", "nu2_1")
    if (length(x$nu1) == 2) rownames(dat) <- c("nu1_1", "nu1_2", "nu2_1", "nu2_2")
    if (length(x$nu1) == 3) rownames(dat) <- c("nu1_1", "nu1_2", "nu1_3", 
                                                  "nu2_1", "nu2_2", "nu2_3")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    cat("\n")
    
    cat("\nRandom effects:                 \n")
    cat("  Formula:", format(as.formula(x$random)), "\n")
    
    if (nrow(x$Sig) == 2) {
      
      var <- all.vars(x$random)[1]
      dat <- matrix(0, nrow = 3, ncol = 4)
      for (i in 1:nrow(x$Sig)) {
        dat[i, ] <- c(x$Sig[i,i], x$seSig[i,i], x$Sig[i,i]/x$seSig[i,i], 2 * pnorm(-abs(x$Sig[i,i]/x$seSig[i,i])))
      }
      dat[3, ] <- c(x$Sig[1,2], x$seSig[1,2], x$Sig[1,2]/x$seSig[1,2], 2 * pnorm(-abs(x$Sig[1,2]/x$seSig[1,2])))
      dat <- as.data.frame(dat)
      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
      name1 <- paste0("(Intercept):", var)
      rownames(dat) <- c("(Intercept)", var, name1)
      dat[, 1:3] <- round(dat[, 1:3], digits+1)
      dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
      print(dat)
      
    } else if (nrow(x$Sig) == 1) {
      
      dat <- matrix(0, nrow = 1, ncol = 4)
      dat[1, ] <- c(x$Sig[1,1], x$seSig[1,1], x$Sig[1,1]/x$seSig[1,1], 2 * pnorm(-abs(x$Sig[1,1]/x$seSig[1,1])))
      
      dat <- as.data.frame(dat)

      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
      rownames(dat) <- "(Intercept)"
      dat[, 1:3] <- round(dat[, 1:3], digits+1)
      dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
      print(dat)
    } else {
      var <- all.vars(x$random)[1:2]
      dat <- matrix(0, nrow = 6, ncol = 4)
      for (i in 1:nrow(x$Sig)) {
        dat[i, ] <- c(x$Sig[i,i], x$seSig[i,i], x$Sig[i,i]/x$seSig[i,i], 2 * pnorm(-abs(x$Sig[i,i]/x$seSig[i,i])))
      }
      dat[4, ] <- c(x$Sig[1,2], x$seSig[1,2], x$Sig[1,2]/x$seSig[1,2], 2 * pnorm(-abs(x$Sig[1,2]/x$seSig[1,2])))
      dat[5, ] <- c(x$Sig[2,3], x$seSig[2,3], x$Sig[2,3]/x$seSig[2,3], 2 * pnorm(-abs(x$Sig[2,3]/x$seSig[2,3])))
      dat[6, ] <- c(x$Sig[1,3], x$seSig[1,3], x$Sig[1,3]/x$seSig[1,3], 2 * pnorm(-abs(x$Sig[1,3]/x$seSig[1,3])))
      dat <- as.data.frame(dat)
      
      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
      rownames(dat) <- c("(Intercept)", var, paste0("(Intercept):", var[1]), 
                         paste0(var[1], ":", var[2]), paste0("(Intercept):", var[2]))                                                            
      dat[, 1:3] <- round(dat[, 1:3], digits+1)
      dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
      print(dat)
      }
  } else {
    cat("Data Summary:\n")
    cat("Number of observations:", nrow(x$ydata), "\n")
    cat("Number of groups:", nrow(x$cdata), "\n\n")
    cat("Proportion of events:", round(x$PropEventType[2, 2]/nrow(x$cdata)*100, 2), "%\n")
    cat("\nNumerical intergration:\n")
    cat(paste("Method:", x$Quad.method, "Guass-Hermite quadrature\n"))
    cat("Number of quadrature points: ", x$quadpoint, "\n")
    cat("\nModel Type: joint modeling of longitudinal continuous and survival data", "\n\n")
    cat("Model summary:\n")
    cat("Longitudinal process: linear mixed effects model\n")
    cat("Event process: Cox proportional hazard model with non-parametric baseline hazard\n\n")
    cat("Loglikelihood: ", x$loglike, "\n\n")
    cat("Fixed effects in the longitudinal submodel: ",
        sprintf(format(paste(deparse(x$LongitudinalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    
    dat <- data.frame(x$beta, x$sebeta, x$beta/x$sebeta, 2 * pnorm(-abs(x$beta/x$sebeta)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    cat("\n")
    dat <- data.frame(x$sigma, x$sesigma, x$sigma/x$sesigma, 2 * pnorm(-abs(x$sigma/x$sesigma)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    rownames(dat) <- "sigma^2"
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
    
    cat("\n Association parameters:                 \n")
    dat <- data.frame(x$nu1, x$senu1, x$nu1/x$senu1, 2 * pnorm(-abs(x$nu1/x$senu1)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    if (length(x$nu1) == 1) rownames(dat) <- c("nu1_1")
    if (length(x$nu1) == 2) rownames(dat) <- c("nu1_1", "nu1_2")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    cat("\n")
    
    cat("\nRandom effects:                 \n")
    cat("  Formula:", format(as.formula(x$random)), "\n")
    
    
    if (nrow(x$Sig) == 2) {
      
      var <- all.vars(x$random)[1]
      dat <- matrix(0, nrow = 3, ncol = 4)
      for (i in 1:nrow(x$Sig)) {
        dat[i, ] <- c(x$Sig[i,i], x$seSig[i,i], x$Sig[i,i]/x$seSig[i,i], 2 * pnorm(-abs(x$Sig[i,i]/x$seSig[i,i])))
      }
      dat[3, ] <- c(x$Sig[1,2], x$seSig[1,2], x$Sig[1,2]/x$seSig[1,2], 2 * pnorm(-abs(x$Sig[1,2]/x$seSig[1,2])))
      dat <- as.data.frame(dat)
      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
      name1 <- paste0("(Intercept):", var)
      rownames(dat) <- c("(Intercept)", var, name1)
      dat[, 1:3] <- round(dat[, 1:3], digits+1)
      dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
      print(dat)
      
    } else if (nrow(x$Sig) == 1) {
      
      dat <- matrix(0, nrow = 1, ncol = 4)
      dat[1, ] <- c(x$Sig[1,1], x$seSig[1,1], x$Sig[1,1]/x$seSig[1,1], 2 * pnorm(-abs(x$Sig[1,1]/x$seSig[1,1])))
      
      dat <- as.data.frame(dat)
      
      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
      rownames(dat) <- "(Intercept)"
      dat[, 1:3] <- round(dat[, 1:3], digits+1)
      dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
      print(dat)
    } else {
      var <- all.vars(x$random)[1:2]
      dat <- matrix(0, nrow = 6, ncol = 4)
      for (i in 1:nrow(x$Sig)) {
        dat[i, ] <- c(x$Sig[i,i], x$seSig[i,i], x$Sig[i,i]/x$seSig[i,i], 2 * pnorm(-abs(x$Sig[i,i]/x$seSig[i,i])))
      }
      dat[4, ] <- c(x$Sig[1,2], x$seSig[1,2], x$Sig[1,2]/x$seSig[1,2], 2 * pnorm(-abs(x$Sig[1,2]/x$seSig[1,2])))
      dat[5, ] <- c(x$Sig[2,3], x$seSig[2,3], x$Sig[2,3]/x$seSig[2,3], 2 * pnorm(-abs(x$Sig[2,3]/x$seSig[2,3])))
      dat[6, ] <- c(x$Sig[1,3], x$seSig[1,3], x$Sig[1,3]/x$seSig[1,3], 2 * pnorm(-abs(x$Sig[1,3]/x$seSig[1,3])))
      dat <- as.data.frame(dat)
      
      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
      rownames(dat) <- c("(Intercept)", var, paste0("(Intercept):", var[1]), 
                         paste0(var[1], ":", var[2]), paste0("(Intercept):", var[2]))                                                            
      dat[, 1:3] <- round(dat[, 1:3], digits+1)
      dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
      print(dat)
    }
    
  }
}
