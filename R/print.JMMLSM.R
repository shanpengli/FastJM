##' Print contents of JMMLSM object.
##' @title Print JMMLSM
##' @param x Object of class 'JMMLSM'.
##' @param digits number of digits of decimal to be printed.
##' @param ... Further arguments passed to or from other methods.
##' @return a summary of data, joint model, log likelihood, and parameter estimates.
##' @seealso \code{\link{JMMLSM}}
##' @author Shanpeng Li
##' @export
##' 
print.JMMLSM <- function(x, digits = 4, ...) {
  if (!inherits(x, "JMMLSM"))
    stop("Use only with 'JMMLSM' objects.\n")
  
  # k = length(x$alpha1)
  # random = all.vars(x$random)  (vector of RE names, no intercept)
  # events = 2 for competing risks; use 1 if single event
  make_assoc_rownames <- function(k, random, events = 2) {
    # Base names: Intercept + random effects
    base <- c("(Intercept)", random)
    # If k exceeds available names, pad with generic RE labels
    if (k > length(base)) {
      base <- c(base, paste0("RE", seq_len(k - length(base))))
    }
    base <- base[seq_len(k)]
    # Attach event suffixes: _1, _2, ...
    unlist(lapply(seq_len(events), function(e) paste0(base, "_", e)))
  }
  
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
    cat("Method: ", x$method, "Guass-Hermite quadrature\n")
    cat("Number of quadrature points: ", x$quadpoint, "\n")
    cat("\nModel Type: joint modeling of longitudinal continuous and competing risks data with the presence of intra-individual variability", "\n\n")
    cat("Model summary:\n")
    cat("Longitudinal process: Mixed effects location scale model\n")
    cat("Event process: cause-specific Cox proportional hazard model with non-parametric baseline hazard\n\n")
    cat("Loglikelihood: ", x$loglike, "\n\n")
    cat("Fixed effects in mean of longitudinal submodel: ",
        sprintf(format(paste(deparse(x$LongitudinalSubmodelmean, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    
    dat <- data.frame(x$beta, x$sebeta, x$beta/x$sebeta, 2 * pnorm(-abs(x$beta/x$sebeta)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)

    cat("\nFixed effects in variance of longitudinal submodel: ",
        sprintf(format(paste(deparse(x$LongitudinalSubmodelvariance, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    
    dat <- data.frame(x$tau, x$setau, x$tau/x$setau, 2 * pnorm(-abs(x$tau/x$setau)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)

    cat("\nSurvival sub-model fixed effects: ",
        sprintf(format(paste(deparse(x$SurvivalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    dat <- data.frame(x$gamma1, x$segamma1, x$gamma1/x$segamma1, 2 * pnorm(-abs(x$gamma1/x$segamma1)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    
    dat <- data.frame(x$gamma2, x$segamma2, x$gamma2/x$segamma2, 2 * pnorm(-abs(x$gamma2/x$segamma2)))
    colnames(dat) <- NULL
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    cat("\nAssociation parameters:                 \n")
    dat <- data.frame(x$alpha1, x$sealpha1, x$alpha1/x$sealpha1, 2 * pnorm(-abs(x$alpha1/x$sealpha1)))
    subdat <- data.frame(x$alpha2, x$sealpha2, x$alpha2/x$sealpha2, 2 * pnorm(-abs(x$alpha2/x$sealpha2)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    colnames(subdat) <- c("Estimate", "SE", "Z value", "p-val")
    dat <- rbind(dat, subdat)
    random <- all.vars(x$random)
    # Usage:
    k <- length(x$alpha1)
    rownames(dat) <- make_assoc_rownames(k, random, events = if (isTRUE(x$CompetingRisk)) 2 else 1)
    
    datnu <- matrix(0, nrow = 2, ncol = 4)
    datnu[1, ] <- c(x$vee1, x$sevee1, x$vee1/x$sevee1, 2 * pnorm(-abs(x$vee1/x$sevee1)))
    datnu[2, ] <- c(x$vee2, x$sevee2, x$vee2/x$sevee2, 2 * pnorm(-abs(x$vee2/x$sevee2)))
    datnu <- as.data.frame(datnu)
    colnames(datnu) <- c("Estimate", "SE", "Z value", "p-val")
    rownames(datnu) <- c("var_(Intercept)_1", "var_(Intercept)_2")
    dat <- rbind(dat, datnu)
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    cat("\n")
    
    cat("\nRandom effects:                 \n")
    cat("  Formula:", format(as.formula(x$random)), "\n")
    
     if (nrow(x$Sig) == 2) {
      
      dat <- matrix(0, nrow = 3, ncol = 4)
      dat[1, ] <- c(x$Sig[1,1], x$seSig[1,1], x$Sig[1,1]/x$seSig[1,1], 2 * pnorm(-abs(x$Sig[1,1]/x$seSig[1,1])))
      dat[2, ] <- c(x$Sig[2,2], x$seSig[2,2], x$Sig[2,2]/x$seSig[2,2], 2 * pnorm(-abs(x$Sig[2,2]/x$seSig[2,2])))
      dat[3, ] <- c(x$Sig[1,2], x$seSig[1,2], x$Sig[1,2]/x$seSig[1,2], 2 * pnorm(-abs(x$Sig[1,2]/x$seSig[1,2])))
      dat <- as.data.frame(dat)
      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
      rownames(dat) <- c("(Intercept)", "var_(Intercept)", "(Intercept):var_(Intercept)")
      dat[, 1:3] <- round(dat[, 1:3], digits+1)
      dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
      print(dat)
    
    } else {
      random.var <- all.vars(x$random)[-length(all.vars(x$random))]
      dat <- matrix(0, nrow = (1+nrow(x$Sig))*nrow(x$Sig)/2, ncol = 4)
      dat <- as.data.frame(dat)
      colnames(x$Sig) <- rownames(x$Sig) <- c("(Intercept)", random.var, "var_(Intercept)")
      for (i in 1:nrow(x$Sig)) {
        dat[i, ] <- c(x$Sig[i,i], x$seSig[i,i], x$Sig[i,i]/x$seSig[i,i], 2 * pnorm(-abs(x$Sig[i,i]/x$seSig[i,i])))
        rownames(dat)[i] <- colnames(x$Sig)[i]
      }
      count <- 1
      for (i in 1:(nrow(x$Sig)-1)) {
        for (j in 0:(nrow(x$Sig)-i-1)) {
          dat[nrow(x$Sig) + count, ] <- c(x$Sig[j+1,i+j+1], x$seSig[j+1,i+j+1], x$Sig[j+1,i+j+1]/x$seSig[j+1,i+j+1], 
                                        2 * pnorm(-abs(x$Sig[j+1,i+j+1]/x$seSig[j+1,i+j+1])))
          rownames(dat)[nrow(x$Sig) + count] <- paste(rownames(x$Sig)[j+1], colnames(x$Sig)[i+j+1], sep = ":")
          count <- count + 1
        }
      }

      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
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
    cat("Method: ", x$method, "Guass-Hermite quadrature\n")
    cat("Number of quadrature points: ", x$quadpoint, "\n")
    cat("\nModel Type: joint modeling of longitudinal continuous and survival data with the presence of intra-individual variability", "\n\n")
    cat("Model summary:\n")
    cat("Longitudinal process: Mixed effects location scale model\n")
    cat("Event process: Cox proportional hazard model with non-parametric baseline hazard\n\n")
    cat("Loglikelihood: ", x$loglike, "\n\n")
    cat("Fixed effects in mean of longitudinal submodel: ",
        sprintf(format(paste(deparse(x$LongitudinalSubmodelmean, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    
    dat <- data.frame(x$beta, x$sebeta, x$beta/x$sebeta, 2 * pnorm(-abs(x$beta/x$sebeta)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    
    cat("\nFixed effects in variance of longitudinal submodel: \n",
        sprintf(format(paste(deparse(x$LongitudinalSubmodelvariance, width.cutoff = 500), collapse=""))), "\n")
    
    dat <- data.frame(x$tau, x$setau, x$tau/x$setau, 2 * pnorm(-abs(x$tau/x$setau)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    
    cat("\nSurvival sub-model fixed effects: ",
        sprintf(format(paste(deparse(x$SurvivalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    dat <- data.frame(x$gamma1, x$segamma1, x$gamma1/x$segamma1, 2 * pnorm(-abs(x$gamma1/x$segamma1)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    
    cat("\n Association parameters:                 \n")
    dat <- data.frame(x$alpha1, x$sealpha1, x$alpha1/x$sealpha1, 2 * pnorm(-abs(x$alpha1/x$sealpha1)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    random <- all.vars(x$random)
    k <- length(x$alpha1)
    rownames(dat) <- make_assoc_rownames(k, random, events = if (isTRUE(x$CompetingRisk)) 2 else 1)
    
    datnu <- matrix(0, nrow = 1, ncol = 4)
    datnu[1, ] <- c(x$vee1, x$sevee1, x$vee1/x$sevee1, 2 * pnorm(-abs(x$vee1/x$sevee1)))
    datnu <- as.data.frame(datnu)
    colnames(datnu) <- c("Estimate", "SE", "Z value", "p-val")
    rownames(datnu) <- c("var_(Intercept)_1")
    dat <- rbind(dat, datnu)
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    cat("\n")
    
    cat("\nRandom effects:                 \n")
    cat("  Formula:", format(as.formula(x$random)), "\n")
    
    
    if (nrow(x$Sig) == 2) {
      
      dat <- matrix(0, nrow = 3, ncol = 4)
      dat[1, ] <- c(x$Sig[1,1], x$seSig[1,1], x$Sig[1,1]/x$seSig[1,1], 2 * pnorm(-abs(x$Sig[1,1]/x$seSig[1,1])))
      dat[2, ] <- c(x$Sig[2,2], x$seSig[2,2], x$Sig[2,2]/x$seSig[2,2], 2 * pnorm(-abs(x$Sig[2,2]/x$seSig[2,2])))
      dat[3, ] <- c(x$Sig[1,2], x$seSig[1,2], x$Sig[1,2]/x$seSig[1,2], 2 * pnorm(-abs(x$Sig[1,2]/x$seSig[1,2])))
      dat <- as.data.frame(dat)
      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
      rownames(dat) <- c("(Intercept)", "var_(Intercept)", "(Intercept):var_(Intercept)")
      dat[, 1:3] <- round(dat[, 1:3], digits+1)
      dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
      print(dat)
      
    } else {
      random.var <- all.vars(x$random)[-length(all.vars(x$random))]
      dat <- matrix(0, nrow = (1+nrow(x$Sig))*nrow(x$Sig)/2, ncol = 4)
      dat <- as.data.frame(dat)
      colnames(x$Sig) <- rownames(x$Sig) <- c("(Intercept)", random.var, "var_(Intercept)")
      for (i in 1:nrow(x$Sig)) {
        dat[i, ] <- c(x$Sig[i,i], x$seSig[i,i], x$Sig[i,i]/x$seSig[i,i], 2 * pnorm(-abs(x$Sig[i,i]/x$seSig[i,i])))
        rownames(dat)[i] <- colnames(x$Sig)[i]
      }
      count <- 1
      for (i in 1:(nrow(x$Sig)-1)) {
        for (j in 0:(nrow(x$Sig)-i-1)) {
          dat[nrow(x$Sig) + count, ] <- c(x$Sig[j+1,i+j+1], x$seSig[j+1,i+j+1], x$Sig[j+1,i+j+1]/x$seSig[j+1,i+j+1], 
                                          2 * pnorm(-abs(x$Sig[j+1,i+j+1]/x$seSig[j+1,i+j+1])))
          rownames(dat)[nrow(x$Sig) + count] <- paste(rownames(x$Sig)[j+1], colnames(x$Sig)[i+j+1], sep = ":")
          count <- count + 1
        }
      }
      
      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
      dat[, 1:3] <- round(dat[, 1:3], digits+1)
      dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
      print(dat)
      
    }
    
  }
}
