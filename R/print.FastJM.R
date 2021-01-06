##' Print contents of FastJM object.
##'
##'
##' @title Print FastJM
##' @param x Object of class 'FastJM'.
##' @param ... Further arguments passed to or from other methods.
##' @seealso \code{\link{jmcs}}
##' @author Shanpeng Li
##' @export
print.FastJM <- function(x, ...) {
  if (!inherits(x, "FastJM"))
    stop("Not a legitimate \"FastJM\" object")

  if (x$type == "jmcs") {
    cat("Call:\n")
    ##need to add function call (Hong)
    cat("Function\n")
    cat("Data Summary:\n")
    cat("Number of observations:", x$SummaryInfo$Numobs, "\n")
    cat("Number of groups:", x$k, "\n\n")
    cat("Proportion of competing risks: \n")
    for (i in 1:2) {
      cat("Risk", i, ":", x$SummaryInfo$PropComp[i+1], "%\n")
    }
    cat("\nNumerical intergration:\n")
    cat("Method: pseudo-adaptive Guass-Hermite quadrature\n")
    cat("Number of quadrature points: ", x$point, "\n")
    cat("\nModel Type: joint modeling of longitudinal continuous and competing risks data", "\n\n")
    cat("Model summary:\n")
    cat("Longitudinal process: linear mixed effects model\n")
    cat("Event process: cause-specific Cox proportional hazard model with non-parametric baseline hazard\n\n")
    cat("Loglikelihood: ", x$loglike, "\n\n")
    cat("Longitudinal sub-model fixed effects: ",
        sprintf(format(paste(deparse(x$SummaryInfo$LongitudinalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n                  Estimate   Std. Error       95% CI                Pr(>|Z|)    \n")
    cat("Longitudinal:                \n")
    cat(" Fixed effects:                 \n")


    for (i in 1:length(x$betas)) {
      #beta = paste0("beta", i)
      beta=names(x$betas)[i]
      uppsd = x$betas[i] + 1.96 * x$se_betas[i]
      lowersd = x$betas[i] - 1.96 * x$se_betas[i]
      zval = (x$betas[i]/x$se_betas[i])
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(beta,width=14,flag="-"), sprintf("% 1.4f", x$betas[i]))
      cat("     ", sprintf("% 1.4f", x$se_betas[i]))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }
    cat(" Random effects:                 \n")
    cat(" ",formatC("sigma^2",width=14,flag="-"), sprintf("% 1.4f", x$sigma2_val))
    cat("     ", sprintf("% 1.4f", x$se_sigma2_val))

    uppsd = x$sigma2_val + 1.96 * x$se_sigma2_val
    lowersd = x$sigma2_val - 1.96 * x$se_sigma2_val
    zval = (x$sigma2_val/x$se_sigma2_val)
    pval = 2 * pnorm(-abs(zval))
    cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
    cat("    ", pval)
    cat("\n")


    cat("\nSurvival sub-model fixed effects: ",
        sprintf(format(paste(deparse(x$SummaryInfo$SurvivalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n                  Estimate   Std. Error       95% CI                Pr(>|Z|)    \n")
    cat("Survival:                \n")
    cat(" Fixed effects:                 \n")
    for (i in 1:dim(x$gamma_matrix)[1]) for (j in 1:dim(x$gamma_matrix)[2]) {
      gamma = paste0("gamma", i, j)
      gamma=colnames(x$gamma_matrix)[j]
      gamma=paste0(gamma,'_', i)
      gammaval = x$gamma_matrix[i, j]
      stdgammaval = x$se_gamma_matrix[i, j]
      uppsd = gammaval + 1.96 * stdgammaval
      lowersd = gammaval - 1.96 * stdgammaval
      zval = (gammaval/stdgammaval)
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(gamma,width=14,flag="-"), sprintf("% 1.4f", gammaval))
      cat("     ", sprintf("% 1.4f", stdgammaval))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }
    cat("\n Association parameters:                 \n")
    for (i in 1:length(x$vee1_estimate)) {
      vee1=paste0("vee_", 1, i)
      uppsd = x$vee1_estimate[i] + 1.96 * x$sd_vee1_estimate[i]
      lowersd = x$vee1_estimate[i] - 1.96 * x$sd_vee1_estimate[i]
      zval = (x$vee1_estimate[i]/x$sd_vee1_estimate[i])
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(vee1,width=14,flag="-"), sprintf("% 1.4f", x$vee1_estimate[i]))
      cat("     ", sprintf("% 1.4f", x$sd_vee1_estimate[i]))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }

    for (i in 1:length(x$vee2_estimate)) {
      vee2=paste0("vee_", 2, i)
      uppsd = x$vee2_estimate[i] + 1.96 * x$sd_vee2_estimate[i]
      lowersd = x$vee2_estimate[i] - 1.96 * x$sd_vee2_estimate[i]
      zval = (x$vee2_estimate[i]/x$sd_vee2_estimate[i])
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(vee2,width=14,flag="-"), sprintf("% 1.4f", x$vee2_estimate[i]))
      cat("     ", sprintf("% 1.4f", x$sd_vee2_estimate[i]))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }

    # cat(' Random effects: \n')

    # restore a matrix from its uppertri
    sd_sigmamatrix = matrix(0, dim(x$sigma_matrix)[1], dim(x$sigma_matrix)[1])

    if (dim(x$sigma_matrix)[1] == 1) {
      sd_sigmamatrix[1, 1] <- x$se_sigma
    } else {
      for (t in 1:dim(x$sigma_matrix)[1]) {
        sd_sigmamatrix[t, t] <- x$se_sigma[t]
      }

      for (q in 2:dim(x$sigma_matrix)[1]) {
        for (t in 1:(dim(x$sigma_matrix)[1] + 1 - q)) {
          sd_sigmamatrix[t, q+t-1] <- x$se_sigma[1+t+(q-1)*(dim(x$sigma_matrix)[1]-1)]
        }
      }
    }

    # print sigmabii
    cat("\nRandom effects:                 \n")
    for (i in 1:(dim(x$sigma_matrix)[1]))
      for (j in 1:(dim(x$sigma_matrix)[2])) {
        if (j < i) {
          next
        } else {
          sigma = paste0("sigma_b", i, j)
          sigmaval = x$sigma_matrix[i, j]
          stdsigmaval = sd_sigmamatrix[i, j]
          uppsd = sigmaval + 1.96 * stdsigmaval
          lowersd = sigmaval - 1.96 * stdsigmaval
          zval = (sigmaval/stdsigmaval)
          pval = 2 * pnorm(-abs(zval))
          cat(" ",formatC(sigma,width=14,flag="-"), sprintf("% 1.4f", sigmaval))
          cat("     ", sprintf("% 1.4f", stdsigmaval))
          cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
          cat("     ", pval)
          cat("\n")
        }


      }

  }

  if (x$type == "jmcsf") {
    cat("Call:\n")
    ##need to add function call (Hong)
    cat("Function\n")
    cat("Data Summary:\n")
    cat("Number of observations:", x$SummaryInfo$Numobs, "\n")
    cat("Number of groups:", x$k, "\n\n")
    cat("Proportion of events:", x$SummaryInfo$PropComp[1+1], "%\n")
    cat("\nNumerical intergration:\n")
    cat("Method: pseudo-adaptive Guass-Hermite quadrature\n")
    cat("Number of quadrature points: ", x$point, "\n")
    cat("\nModel Type: joint modeling of longitudinal continuous and survival data with single failure type", "\n\n")
    cat("Model summary:\n")
    cat("Longitudinal process: linear mixed effects model\n")
    cat("Event process: cause-specific Cox proportional hazard model with non-parametric baseline hazard\n\n")
    cat("Loglikelihood: ", x$loglike, "\n\n")
    cat("Longitudinal sub-model fixed effects: ",
        sprintf(format(paste(deparse(x$SummaryInfo$LongitudinalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n                  Estimate   Std. Error       95% CI                Pr(>|Z|)    \n")
    cat("Longitudinal:                \n")
    cat(" Fixed effects:                 \n")

    for (i in 1:length(x$betas)) {
      #beta = paste0("beta", i)
      beta=names(x$betas)[i]
      uppsd = x$betas[i] + 1.96 * x$se_betas[i]
      lowersd = x$betas[i] - 1.96 * x$se_betas[i]
      zval = (x$betas[i]/x$se_betas[i])
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(beta,width=14,flag="-"), sprintf("% 1.4f", x$betas[i]))
      cat("     ", sprintf("% 1.4f", x$se_betas[i]))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }
    cat(" Random effects:                 \n")
    cat(" ",formatC("sigma^2",width=14,flag="-"), sprintf("% 1.4f", x$sigma2_val))
    cat("     ", sprintf("% 1.4f", x$se_sigma2_val))

    uppsd = x$sigma2_val + 1.96 * x$se_sigma2_val
    lowersd = x$sigma2_val - 1.96 * x$se_sigma2_val
    zval = (x$sigma2_val/x$se_sigma2_val)
    pval = 2 * pnorm(-abs(zval))
    cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
    cat("    ", pval)
    cat("\n")


    cat("\nSurvival sub-model fixed effects: ",
        sprintf(format(paste(deparse(x$SummaryInfo$SurvivalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n                  Estimate   Std. Error       95% CI                Pr(>|Z|)    \n")
    cat("Survival:                \n")
    cat(" Fixed effects:                 \n")
    for (i in 1:dim(x$gamma_matrix)[1]) for (j in 1:dim(x$gamma_matrix)[2]) {
      gamma = paste0("gamma", i, j)
      gamma=colnames(x$gamma_matrix)[j]
      gamma=paste0(gamma,'_', i)
      gammaval = x$gamma_matrix[i, j]
      stdgammaval = x$se_gamma_matrix[i, j]
      uppsd = gammaval + 1.96 * stdgammaval
      lowersd = gammaval - 1.96 * stdgammaval
      zval = (gammaval/stdgammaval)
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(gamma,width=14,flag="-"), sprintf("% 1.4f", gammaval))
      cat("     ", sprintf("% 1.4f", stdgammaval))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }
    cat("\n Association parameters:                 \n")
    for (i in 1:length(x$vee1_estimate)) {
      vee1=paste0("vee_", 1, i)
      uppsd = x$vee1_estimate[i] + 1.96 * x$sd_vee1_estimate[i]
      lowersd = x$vee1_estimate[i] - 1.96 * x$sd_vee1_estimate[i]
      zval = (x$vee1_estimate[i]/x$sd_vee1_estimate[i])
      pval = 2 * pnorm(-abs(zval))
      cat(" ",formatC(vee1,width=14,flag="-"), sprintf("% 1.4f", x$vee1_estimate[i]))
      cat("     ", sprintf("% 1.4f", x$sd_vee1_estimate[i]))
      cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
      cat("     ", pval)
      cat("\n")
    }

    # cat(' Random effects: \n')

    # restore a matrix from its uppertri
    sd_sigmamatrix = matrix(0, dim(x$sigma_matrix)[1], dim(x$sigma_matrix)[1])

    if (dim(x$sigma_matrix)[1] == 1) {
      sd_sigmamatrix[1, 1] <- x$se_sigma
    } else {
      for (t in 1:dim(x$sigma_matrix)[1]) {
        sd_sigmamatrix[t, t] <- x$se_sigma[t]
      }

      for (q in 2:dim(x$sigma_matrix)[1]) {
        for (t in 1:(dim(x$sigma_matrix)[1] + 1 - q)) {
          sd_sigmamatrix[t, q+t-1] <- x$se_sigma[1+t+(q-1)*(dim(x$sigma_matrix)[1]-1)]
        }
      }
    }

    # print sigmabii
    cat("\nRandom effects:                 \n")
    for (i in 1:(dim(x$sigma_matrix)[1]))
      for (j in 1:(dim(x$sigma_matrix)[2])) {
        if (j < i) {
          next
        } else {
          sigma = paste0("sigma_b", i, j)
          sigmaval = x$sigma_matrix[i, j]
          stdsigmaval = sd_sigmamatrix[i, j]
          uppsd = sigmaval + 1.96 * stdsigmaval
          lowersd = sigmaval - 1.96 * stdsigmaval
          zval = (sigmaval/stdsigmaval)
          pval = 2 * pnorm(-abs(zval))
          cat(" ",formatC(sigma,width=14,flag="-"), sprintf("% 1.4f", sigmaval))
          cat("     ", sprintf("% 1.4f", stdsigmaval))
          cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
          cat("     ", pval)
          cat("\n")
        }


      }

  }
}
