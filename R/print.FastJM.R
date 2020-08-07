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
    cat("Model Type:              ", "Joint modeling with competing risks data", "\n\n")
    cat("                  Estimate   Std. Error       95% CI                Pr(>|Z|)    \n")
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

    cat(" ",formatC("sigma^2",width=14,flag="-"), sprintf("% 1.4f", x$sigma2_val))
    cat("     ", sprintf("% 1.4f", x$se_sigma2_val))

    uppsd = x$sigma2_val + 1.96 * x$se_sigma2_val
    lowersd = x$sigma2_val - 1.96 * x$se_sigma2_val
    zval = (x$sigma2_val/x$se_sigma2_val)
    pval = 2 * pnorm(-abs(zval))
    cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
    cat("    ", pval)
    cat("\n")


    # cat(' Estimate Std. Error 95%CI Pr(>|Z|) \n')
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
    sd_sigmamatrix[lower.tri(sd_sigmamatrix, diag = TRUE)] <- x$se_sigma
    sd_sigmamatrix <- t(sd_sigmamatrix)

    # print sigmabii
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
    cat("Model Type:           ", "Joint modeling with single failure survival data", "\n\n")
    cat("                  Estimate   Std. Error       95% CI                Pr(>|Z|)    \n")
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

    cat(" ",formatC("sigma^2",width=14,flag="-"), sprintf("% 1.4f", x$sigma2_val))
    cat("     ", sprintf("% 1.4f", x$se_sigma2_val))

    uppsd = x$sigma2_val + 1.96 * x$se_sigma2_val
    lowersd = x$sigma2_val - 1.96 * x$se_sigma2_val
    zval = (x$sigma2_val/x$se_sigma2_val)
    pval = 2 * pnorm(-abs(zval))
    cat("     ", sprintf("(% 1.4f,% 1.4f)", lowersd, uppsd))
    cat("    ", pval)
    cat("\n")


    # cat(' Estimate Std. Error 95%CI Pr(>|Z|) \n')
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
    sd_sigmamatrix[lower.tri(sd_sigmamatrix, diag = TRUE)] <- x$se_sigma
    sd_sigmamatrix <- t(sd_sigmamatrix)

    # print sigmabii
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
