##' Joint modelling for longitutal and censored data with competing risks
##' @title Extraction of standard error and 95\% confidence interval of longitudinal/survival sub-model fixed effects
##' @param object  The FastJM object returned by either jmo or jmc function.
##' @param coeff The coefficients returned by the FastJM object. Results of longitudinal/survival sub-model will be printed out when typing "longitudinal" or "survival".
##' @param digits The number of digits to print out.
##' @param ... further arguments passed to or from other methods.
##' @seealso \code{\link{jmc}}
##' @return Return standard errors of parameters with variable names
##' @references
##' \itemize{
##' \item Elashoff, Robert M., Gang Li, and Ning Li. "A joint model for longitudinal measurements and survival data in the presence of multiple failure types." Biometrics 64.3 (2008): 762-771.
##' }
##' @export
summary.FastJM <-
  function (object, coeff=c("longitudinal", "survival"), digits=4, ...) {
    if (!inherits(object, "FastJM"))
      stop("Use only with 'FastJM' objects.\n")
    if (object$type == "jmcs") {
      if (coeff == "longitudinal") {
        ##Estimates of betas
        Estimate <- object$betas
        SE <- object$se_betas
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
      } else if (coeff == "survival") {
        p1a <- length(object$vee1_estimate)
        Estimate <- t(object$gamma_matrix)
        Estimate <- reshape2::melt(Estimate)
        SE <- t(object$se_gamma_matrix)
        SE <- reshape2::melt(SE)
        LowerLimit <- Estimate[, 3] - 1.96 * SE[, 3]
        UpperLimit <- Estimate[, 3] + 1.96 * SE[, 3]
        zval = (Estimate[, 3]/SE[, 3])
        pval = 2 * pnorm(-abs(zval))
        out <- data.frame(Estimate, exp(Estimate[, 3]), SE[, 3], LowerLimit, UpperLimit, pval)
        out[, 1] <- paste(out[, 1], out[, 2], sep = "_")
        out <- out[, -2]
        names(out) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", "p-values")

        ##vee
        Estimate <- c(object$vee1_estimate, object$vee2_estimate)
        SE <- c(object$sd_vee1_estimate, object$sd_vee2_estimate)
        LowerLimit <- Estimate - 1.96 * SE
        UpperLimit <- Estimate + 1.96 * SE
        zval = Estimate/SE
        pval = 2 * pnorm(-abs(zval))
        outvee <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, pval)
        Survival <- vector()
        for (i in 1:p1a) Survival[i] <- paste0("vee_", 1, i)
        for (i in 1:p1a) Survival[i+p1a] <- paste0("vee_", 2, i)
        outvee <- cbind(Survival, outvee)
        ##round p-val (TBD Hong Wang)
        names(outvee) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", "p-values")
        out <- rbind(out, outvee)
        rownames(out) <- NULL
        out[, 2:ncol(out)] <- round(out[, 2:ncol(out)], digits = digits)
        out[, ncol(out)] <- format(out[, ncol(out)], scientific = FALSE)

        return(out)
      }  else {
        stop("Unexpected arguments! Must choose one of the following options: longitudinal, survival")
      }
    } else if (object$type == "jmcsf") {
      if (coeff == "longitudinal") {
        ##Estimates of betas
        Estimate <- object$betas
        SE <- object$se_betas
        LowerLimit <- Estimate - 1.96 * SE
        UpperLimit <- Estimate + 1.96 * SE
        zval = (Estimate/SE)
        pval = 2 * pnorm(-abs(zval))
        out <- data.frame(Estimate, SE, LowerLimit, UpperLimit, pval)
        out <- cbind(rownames(out), out)
        names(out) <- c("Longitudinal", "coef", "SE", "95%Lower", "95%Upper", "p-values")
        rownames(out) <- NULL
        out[, 2:ncol(out)] <- round(out[, 2:ncol(out)], digits = digits)
        out[, ncol(out)] <- format(out[, ncol(out)], scientific = FALSE)

        return(out)
      } else if (coeff == "survival") {
        p1a <- length(object$vee1_estimate)
        Estimate <- t(object$gamma_matrix)
        SE <- t(object$se_gamma_matrix)
        LowerLimit <- Estimate[, 1] - 1.96 * SE[, 1]
        UpperLimit <- Estimate[, 1] + 1.96 * SE[, 1]
        zval = (Estimate[, 1]/SE[, 1])
        pval = 2 * pnorm(-abs(zval))
        out <- data.frame(Estimate[, 1], exp(Estimate[, 1]), SE[, 1], LowerLimit, UpperLimit, pval)
        out <- cbind(rownames(out), out)
        rownames(out) <- NULL
        names(out) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", "p-values")

        ##vee
        Estimate <- object$vee1_estimate
        SE <- object$sd_vee1_estimate
        LowerLimit <- Estimate - 1.96 * SE
        UpperLimit <- Estimate + 1.96 * SE
        zval = Estimate/SE
        pval = 2 * pnorm(-abs(zval))
        outvee <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, pval)
        Survival <- vector()
        for (i in 1:p1a) Survival[i] <- paste0("vee_", 1, i)
        outvee <- cbind(Survival, outvee)
        ##round p-val (TBD Hong Wang)
        names(outvee) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", "p-values")
        out <- rbind(out, outvee)
        rownames(out) <- NULL
        out[, 2:ncol(out)] <- round(out[, 2:ncol(out)], digits = digits)
        out[, ncol(out)] <- format(out[, ncol(out)], scientific = FALSE)

        return(out)
      }  else {
        stop("Unexpected arguments! Must choose one of the following options: longitudinal, survival")
      }
    }
  }
