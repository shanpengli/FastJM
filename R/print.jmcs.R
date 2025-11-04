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
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    cat("\n")
    
    cat("\nRandom effects:                 \n")
    cat("  Formula:", format(as.formula(x$random)), "\n")
    
    stopifnot(is.matrix(x$Sig), nrow(x$Sig) == ncol(x$Sig))
    p <- nrow(x$Sig)
    stopifnot(is.matrix(x$seSig), all(dim(x$seSig) == c(p, p)))
    
    # Build effect names: "(Intercept)", then variables from x$random (up to p-1)
    vars <- all.vars(x$random)
    eff_names <- c("(Intercept)", head(vars, max(0, p - 1)))
    if (length(eff_names) < p) {
      # pad safely if fewer names available
      eff_names <- c(eff_names, paste0("RE", seq_len(p - length(eff_names)) + (length(eff_names) > 0)))
    }
    eff_names <- eff_names[seq_len(p)]
    
    # Row index sets: first variances (i,i), then covariances (i<j)
    var_idx <- cbind(seq_len(p), seq_len(p))
    cov_idx <- t(combn(p, 2))  # each row = (i, j) with i<j
    
    idx <- rbind(var_idx, cov_idx)
    n_rows <- nrow(idx)
    
    # Compute estimates, SE, Z, p
    dat <- matrix(NA_real_, nrow = n_rows, ncol = 4)
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    
    for (r in seq_len(n_rows)) {
      i <- idx[r, 1]; j <- idx[r, 2]
      est <- x$Sig[i, j]
      se  <- x$seSig[i, j]
      z   <- if (is.finite(se) && se > 0) est / se else NA_real_
      pval <- if (is.finite(z)) 2 * pnorm(-abs(z)) else NA_real_
      dat[r, ] <- c(est, se, z, pval)
    }
    
    # Row names: variances get eff_names; covariances get "name_i:name_j"
    rownms <- c(eff_names, paste0(eff_names[cov_idx[,1]], ":", eff_names[cov_idx[,2]]))
    dat <- as.data.frame(dat, row.names = rownms)
    
    # Formatting
    dat[, 1:3] <- round(dat[, 1:3], digits + 1)
    dat[, 4]   <- sprintf(paste0("%.", digits, "f"), dat[, 4])
    
    print(dat)
    
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
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    cat("\n")
    
    cat("\nRandom effects:                 \n")
    cat("  Formula:", format(as.formula(x$random)), "\n")
    
    stopifnot(is.matrix(x$Sig), nrow(x$Sig) == ncol(x$Sig))
    p <- nrow(x$Sig)
    stopifnot(is.matrix(x$seSig), all(dim(x$seSig) == c(p, p)))
    
    # Build effect names: "(Intercept)", then variables from x$random (up to p-1)
    vars <- all.vars(x$random)
    eff_names <- c("(Intercept)", head(vars, max(0, p - 1)))
    if (length(eff_names) < p) {
      # pad safely if fewer names available
      eff_names <- c(eff_names, paste0("RE", seq_len(p - length(eff_names)) + (length(eff_names) > 0)))
    }
    eff_names <- eff_names[seq_len(p)]
    
    # Row index sets: first variances (i,i), then covariances (i<j)
    var_idx <- cbind(seq_len(p), seq_len(p))
    cov_idx <- t(combn(p, 2))  # each row = (i, j) with i<j
    
    idx <- rbind(var_idx, cov_idx)
    n_rows <- nrow(idx)
    
    # Compute estimates, SE, Z, p
    dat <- matrix(NA_real_, nrow = n_rows, ncol = 4)
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    
    for (r in seq_len(n_rows)) {
      i <- idx[r, 1]; j <- idx[r, 2]
      est <- x$Sig[i, j]
      se  <- x$seSig[i, j]
      z   <- if (is.finite(se) && se > 0) est / se else NA_real_
      pval <- if (is.finite(z)) 2 * pnorm(-abs(z)) else NA_real_
      dat[r, ] <- c(est, se, z, pval)
    }
    
    # Row names: variances get eff_names; covariances get "name_i:name_j"
    rownms <- c(eff_names, paste0(eff_names[cov_idx[,1]], ":", eff_names[cov_idx[,2]]))
    dat <- as.data.frame(dat, row.names = rownms)
    
    # Formatting
    dat[, 1:3] <- round(dat[, 1:3], digits + 1)
    dat[, 4]   <- sprintf(paste0("%.", digits, "f"), dat[, 4])
    
    print(dat)
    
  }
}
