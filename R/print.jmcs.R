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
    dat[, 1:3] <- round(dat[, 1:3], digits)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    sigma2 <- as.numeric(x$sigma)
    residual_tab <- data.frame(
      Variance = sigma2,
      StdDev = sqrt(sigma2),
      row.names = "Residual"
    )
    residual_tab[] <- lapply(
      residual_tab,
      round,
      digits = digits
    )
    cat("\nResidual error:\n")
    print(residual_tab)
    
    cat("\nFixed effects in the survival sub-model: ",
        sprintf(format(paste(deparse(x$SurvivalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    dat <- data.frame(x$gamma1, x$segamma1, x$gamma1/x$segamma1, 2 * pnorm(-abs(x$gamma1/x$segamma1)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    
    dat2 <- data.frame(x$gamma2, x$segamma2, x$gamma2/x$segamma2, 2 * pnorm(-abs(x$gamma2/x$segamma2)))
    colnames(dat2) <- c("Estimate", "SE", "Z value", "p-val")
    dat <- rbind(dat, dat2)
    dat[, 1:3] <- round(dat[, 1:3], digits)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    cat("\nAssociation parameters:                 \n")
    dat <- data.frame(x$nu1, x$senu1, x$nu1/x$senu1, 2 * pnorm(-abs(x$nu1/x$senu1)))
    subdat <- data.frame(x$nu2, x$senu2, x$nu2/x$senu2, 2 * pnorm(-abs(x$nu2/x$senu2)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    colnames(subdat) <- c("Estimate", "SE", "Z value", "p-val")
    dat <- rbind(dat, subdat)
    random <- all.vars(x$random)
    
    q <- length(x$nu1)
    stopifnot(q >= 1)
    # Infer number of outcomes (blocks) from dat
    if (nrow(dat) %% q != 0) stop("nrow(dat) must be a multiple of length(x$nu1).")
    G <- nrow(dat) %/% q
    # We need q - 1 names beyond the intercept
    need <- max(0, q - 1)
    if (length(random) < need) {
      stop("`random` must have at least length(x$nu1) - 1 names.")
    }
    base <- c("(Intercept)", random[seq_len(need)])
    rownames(dat) <- unlist(lapply(seq_len(G), function(k) paste0(base, "_", k)))
    dat[, 1:3] <- round(dat[, 1:3], digits)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    cat("\nRandom effects:                 \n")
    cat("  Formula:", format(as.formula(x$random)), "\n")
    vars <- all.vars(x$random)
    p <- nrow(x$Sig)
    Sig <- as.matrix(x$Sig)
    Sig <- (Sig + t(Sig)) / 2
    
    re_sd  <- sqrt(pmax(diag(Sig), 0))
    re_cor <- cov2cor(Sig)
    
    vc <- matrix("", nrow = p, ncol = p)
    
    eff_names <- c("(Intercept)", head(vars, max(0, p - 1)))
    rownames(vc) <- eff_names
    
    corr_names <- eff_names[seq_len(p - 1)]
    corr_names[corr_names == "(Intercept)"] <- "(Intr)"
    
    colnames(vc) <- c(
      "StdDev",
      if (p > 1) corr_names else NULL
    )
    
    vc[, 1] <- formatC(
      re_sd,
      format = "f",
      digits = digits
    )
    
    if (p > 1) {
      for (i in 2:p) {
        vc[i, 2:i] <- formatC(
          re_cor[i, seq_len(i - 1)],
          format = "f",
          digits = digits
        )
      }
    }
    print(vc, quote = FALSE, right = TRUE)
    
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
    dat[, 1:3] <- round(dat[, 1:3], digits)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    sigma2 <- as.numeric(x$sigma)
    residual_tab <- data.frame(
      Variance = sigma2,
      StdDev = sqrt(sigma2),
      row.names = "Residual"
    )
    residual_tab[] <- lapply(
      residual_tab,
      round,
      digits = digits
    )
    cat("\nResidual error:\n")
    print(residual_tab)
    cat("\nFixed effects in the survival sub-model: ",
        sprintf(format(paste(deparse(x$SurvivalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    dat <- data.frame(x$gamma1, x$segamma1, x$gamma1/x$segamma1, 2 * pnorm(-abs(x$gamma1/x$segamma1)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    
    cat("\n Association parameters:                 \n")
    dat <- data.frame(x$nu1, x$senu1, x$nu1/x$senu1, 2 * pnorm(-abs(x$nu1/x$senu1)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    random <- all.vars(x$random)
    q <- length(x$nu1)
    stopifnot(q >= 1)
    # Infer number of outcomes (blocks) from dat
    if (nrow(dat) %% q != 0) stop("nrow(dat) must be a multiple of length(x$nu1).")
    G <- nrow(dat) %/% q
    # We need q - 1 names beyond the intercept
    need <- max(0, q - 1)
    if (length(random) < need) {
      stop("`random` must have at least length(x$nu1) - 1 names.")
    }
    base <- c("(Intercept)", random[seq_len(need)])
    rownames(dat) <- unlist(lapply(seq_len(G), function(k) paste0(base, "_", k)))
    dat[, 1:3] <- round(dat[, 1:3], digits)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    cat("\n")
    
    cat("\nRandom effects:                 \n")
    cat("  Formula:", format(as.formula(x$random)), "\n")
    vars <- all.vars(x$random)
    p <- nrow(x$Sig)
    Sig <- as.matrix(x$Sig)
    Sig <- (Sig + t(Sig)) / 2
    
    re_sd  <- sqrt(pmax(diag(Sig), 0))
    re_cor <- cov2cor(Sig)
    
    vc <- matrix("", nrow = p, ncol = p)
    eff_names <- c("(Intercept)", head(vars, max(0, p - 1)))
    rownames(vc) <- eff_names
    
    corr_names <- eff_names[seq_len(p - 1)]
    corr_names[corr_names == "(Intercept)"] <- "(Intr)"
    
    colnames(vc) <- c(
      "StdDev",
      if (p > 1) corr_names else NULL
    )
    
    vc[, 1] <- formatC(
      re_sd,
      format = "f",
      digits = digits
    )
    
    if (p > 1) {
      for (i in 2:p) {
        vc[i, 2:i] <- formatC(
          re_cor[i, seq_len(i - 1)],
          format = "f",
          digits = digits
        )
      }
    }
    print(vc, quote = FALSE, right = TRUE)
    
  }
}
