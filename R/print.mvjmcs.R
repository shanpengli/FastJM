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
    if (x$latAsso == "sre") {
      cat("Number of observations:", nrow(x$ydata), "\n")
      cat("Number of groups:", nrow(x$cdata), "\n\n")
    } else {
      cdata <- x$cdata
      ydata <- x$ydata
      s <- x$s
      time <- all.vars(x$SurvivalSubmodel)[1]
      random <- x$random
      if (is.list(random)) {
        random.form <- all.vars(random[[1]])
        ID <- random.form[length(random.form)]
      } else {
        random.form <- all.vars(random)
        ID <- random.form[length(random.form)]
      }
      subcdata <- cdata %>%
        dplyr::filter(.data[[time]] > s)
      cID <- subcdata[, ID]
      subydata <- ydata %>%
        dplyr::filter(.data[[ID]] %in% cID)
      cat("Number of observations:", nrow(subydata), "\n")
      cat("Number of groups:", nrow(subcdata), "\n\n")
    }
    cat("Proportion of competing risks: \n")
    if (x$latAsso == "sre") {
      for (i in 1:2) {
        cat("Risk", i, ":", round(x$PropEventType[i+1, 2]/nrow(x$cdata)*100, 2), "%\n")
      }
    } else {
      status <- all.vars(x$SurvivalSubmodel)[2]
      PropComp <- as.data.frame(table(subcdata[, status]))
      for (i in 1:2) {
        cat("Risk", i, ":", round(PropComp[i+1, 2]/nrow(subcdata)*100, 2), "%\n")
      }
    }
    cat("\nModel Type: joint modeling of multivariate longitudinal continuous and competing risks data", "\n\n")
    cat("Model summary:\n")
    if (x$latAsso != "sre") {
      cat(sprintf("Landmark analysis: Yes (s = %s)\n", s))
      if (x$latAsso == "presentlp") cat("Latent association: current value of the latent process\n")
      if (x$latAsso == "present") cat("Latent association: current value\n")
    }
    if (!is.null(x$runtime)) {
      cat("Runtime:", format_runtime(x$runtime), "\n")
    }
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
    dat[, 1:3] <- round(dat[, 1:3], digits)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    cat("\n")
    
    # ~~~~~~~~~~~~~~~~~
    # print sigma^2 est
    # ~~~~~~~~~~~~~~~~~
    sigma2 <- as.numeric(x$sigma)
    residual_sd <- sqrt(sigma2)

    residual_tab <- data.frame(
      Variance = sigma2,
      StdDev = residual_sd,
      row.names = paste0("sigma_bio", seq_len(numBio)),
      check.names = FALSE
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
    
    # ~~~~~~~~~~~~~~~~~
    # print gamma est
    # ~~~~~~~~~~~~~~~~~
    gamma1 <- x$gamma1
    segamma1 <- x$segamma1
    gamma2 <- x$gamma2
    segamma2 <- x$segamma2
    
    dat <- data.frame(
      gamma1,
      segamma1,
      gamma1 / segamma1,
      2 * pnorm(-abs(gamma1 / segamma1))
    )
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    
    dat2 <- data.frame(
      gamma2,
      segamma2,
      gamma2 / segamma2,
      2 * pnorm(-abs(gamma2 / segamma2))
    )
    colnames(dat2) <- c("Estimate", "SE", "Z value", "p-val")
    
    dat <- rbind(dat, dat2)
    dat[, 1:3] <- round(dat[, 1:3], digits)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    # ~~~~~~~~~~~~~~~~~
    # print association est
    # ~~~~~~~~~~~~~~~~~
    
    cat("\nAssociation parameters:                 \n")
    dat <- data.frame(x$alpha1, x$sealpha1, x$alpha1/x$sealpha1, 2 * pnorm(-abs(x$alpha1/x$sealpha1)))
    subdat <- data.frame(x$alpha2, x$sealpha2, x$alpha2/x$sealpha2, 2 * pnorm(-abs(x$alpha2/x$sealpha2)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    colnames(subdat) <- c("Estimate", "SE", "Z value", "p-val")
    dat <- rbind(dat, subdat)
    tempName <- c()
    ind = 1
    
    for(g in 1:numBio){
      if (x$latAsso == "sre") {
        
        pRE <- length(all.vars(x$random[[g]]))
        temp <- all.vars(x$random[[g]])
        
        tempName[ind:(ind+pRE-1)] <-
          paste0(c("(Intercept)", temp[-length(temp)]), "_1bio", g)
        
      } else if (x$latAsso %in% c("present", "presentlp")) {
        
        pRE <- 1
        tempName[ind] <- paste0("alpha1_bio", g)
        
      } else {
        stop("Unknown latent association structure.")
      }
      
      ind <- ind + pRE
    }
    
    for(g in 1:numBio){
      
      if (x$latAsso == "sre") {
        
        pRE <- length(all.vars(x$random[[g]]))
        temp <- all.vars(x$random[[g]])
        
        tempName[ind:(ind+pRE-1)] <-
          paste0(c("(Intercept)", temp[-length(temp)]), "_2bio", g)
        
      } else if (x$latAsso %in% c("present", "presentlp")) {
        
        pRE <- 1
        tempName[ind] <- paste0("alpha2_bio", g)
        
      } else {
        stop("Unknown latent association structure.")
      }
      
      ind <- ind + pRE
      
      
    }
    
    rownames(dat) <- tempName
    
    
    
    dat[, 1:3] <- round(dat[, 1:3], digits)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    cat("\n")
    
    cat("\nRandom effects:                 \n")
    for(g in 1:numBio){
      cat("  bio",g, ": ", format(as.formula(x$random[[g]])), "\n")
    }
    pREtotal <- nrow(x$Sig)
    
    Sig <- as.matrix(x$Sig)
    
    # Protect against negligible numerical asymmetry
    Sig <- (Sig + t(Sig)) / 2
    
    # ---------------------------------------------------------
    # Construct random-effect names for all biomarkers
    # ---------------------------------------------------------
    
    re_names <- character(0)
    
    for (g in seq_len(numBio)) {
      
      random_vars <- all.vars(x$random[[g]])
      
      if (length(random_vars) > 1) {
        random_vars <- random_vars[-length(random_vars)]
        biomarker_names <- c(
          paste0("Intercept", g),
          paste0(random_vars, g)
        )
      } else {
        biomarker_names <- paste0("Intercept", g)
      }
      re_names <- c(re_names, biomarker_names)
    }
    
    if (length(re_names) != pREtotal) {
      stop(
        "The number of constructed random-effect names (",
        length(re_names),
        ") does not match the dimension of x$Sig (",
        pREtotal,
        ")."
      )
    }
    
    # ---------------------------------------------------------
    # Calculate standard deviations and correlations
    # ---------------------------------------------------------
    
    if (any(diag(Sig) < 0)) {
      warning(
        "Negative diagonal elements were found in x$Sig; ",
        "the corresponding standard deviations are set to zero."
      )
    }
    
    re_sd <- sqrt(pmax(diag(Sig), 0))
    
    # Calculate correlations safely, including zero-variance cases
    denom <- outer(re_sd, re_sd)
    re_cor <- Sig / denom
    
    re_cor[!is.finite(re_cor)] <- NA_real_
    diag(re_cor) <- ifelse(re_sd > 0, 1, NA_real_)
    
    vc <- matrix(
      "",
      nrow = pREtotal,
      ncol = pREtotal
    )
    
    rownames(vc) <- re_names
    
    corr_names <- if (pREtotal > 1) {
      re_names[seq_len(pREtotal - 1)]
    } else {
      character(0)
    }
    
    # Shorter correlation-column labels
    corr_names <- sub("^Intercept", "Intr", corr_names)
    
    colnames(vc) <- c(
      "StdDev",
      corr_names
    )
    
    # Standard deviations
    vc[, 1] <- formatC(
      re_sd,
      format = "f",
      digits = digits
    )
    
    # Lower-triangular correlations
    if (pREtotal > 1) {
      for (i in 2:pREtotal) {
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
    if (x$latAsso == "sre") {
      cat("Number of observations:", nrow(x$ydata), "\n")
      cat("Number of groups:", nrow(x$cdata), "\n\n")
    } else {
      cdata <- x$cdata
      ydata <- x$ydata
      s <- x$s
      time <- all.vars(x$SurvivalSubmodel)[1]
      random <- x$random
      if (is.list(random)) {
        random.form <- all.vars(random[[1]])
        ID <- random.form[length(random.form)]
      } else {
        random.form <- all.vars(random)
        ID <- random.form[length(random.form)]
      }
      subcdata <- cdata %>%
        dplyr::filter(.data[[time]] > s)
      cID <- subcdata[, ID]
      subydata <- ydata %>%
        dplyr::filter(.data[[ID]] %in% cID)
      cat("Number of observations:", nrow(subydata), "\n")
      cat("Number of groups:", nrow(subcdata), "\n\n")
    }
    if (x$latAsso == "sre") {
      cat("Proportion of events:", round(x$PropEventType[2, 2]/nrow(x$cdata)*100, 2), "%\n")
    } else {
      status <- all.vars(x$SurvivalSubmodel)[2]
      PropComp <- as.data.frame(table(subcdata[, status]))
      cat("Proportion of events:", round(PropComp[2, 2]/nrow(subcdata)*100, 2), "%\n")
    }
    cat("\nModel Type: joint modeling of multivariate longitudinal continuous and survival data", "\n\n")
    cat("Model summary:\n")
    if (x$latAsso != "sre") {
      cat(sprintf("Landmark analysis: Yes (s = %s)\n", s))
      if (x$latAsso == "presentlp") cat("Latent association: current value of the latent process\n")
      if (x$latAsso == "present") cat("Latent association: current value\n")
    }
    if (!is.null(x$runtime)) {
      cat("Runtime:", format_runtime(x$runtime), "\n")
    }
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
    dat[, 1:3] <- round(dat[, 1:3], digits)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    cat("\n")
    # ~~~~~~~~~~~~~~~~~
    # print sigma^2 est
    # ~~~~~~~~~~~~~~~~~
    sigma2 <- as.numeric(x$sigma)
    residual_sd <- sqrt(sigma2)
    
    residual_tab <- data.frame(
      Variance = sigma2,
      StdDev = residual_sd,
      row.names = paste0("sigma_bio", seq_len(numBio)),
      check.names = FALSE
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
    
    # ~~~~~~~~~~~~~~~~~
    # print gamma est
    # ~~~~~~~~~~~~~~~~~
    gamma1 <- x$gamma1
    segamma1 <- x$segamma1
    dat <- data.frame(
      gamma1,
      segamma1,
      gamma1 / segamma1,
      2 * pnorm(-abs(gamma1 / segamma1))
    )
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    
    dat[, 1:3] <- round(dat[, 1:3], digits)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
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
      if (x$latAsso == "sre") {
        
        pRE <- length(all.vars(x$random[[g]]))
        temp <- all.vars(x$random[[g]])
        
        if (pRE == 1){
          
          tempName[ind] <- paste0("(Intercept)_1bio", g)
          
        } else {
          
          tempName[ind:(ind+pRE-1)] <-
            paste0(
              c("(Intercept)", temp[-length(temp)]),
              "_1bio", g
            )
        }
        
      } else if (x$latAsso %in% c("pv", "pvlp", "present", "presentlp")) {
        
        pRE <- 1
        tempName[ind] <- paste0("alpha1_bio", g)
        
      } else {
        
        stop("Unknown latent association structure.")
        
      }
      
      ind <- ind + pRE
    }
    
    rownames(dat) <- tempName
    
    dat[, 1:3] <- round(dat[, 1:3], digits)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    cat("\n")
    
    cat("\nRandom effects:                 \n")
    for(g in 1:numBio){
      cat("  bio",g, ": ", format(as.formula(x$random[[g]])), "\n")
    }
    pREtotal <- nrow(x$Sig)
    
    Sig <- as.matrix(x$Sig)
    
    # Protect against negligible numerical asymmetry
    Sig <- (Sig + t(Sig)) / 2
    
    # ---------------------------------------------------------
    # Construct random-effect names for all biomarkers
    # ---------------------------------------------------------
    
    re_names <- character(0)
    
    for (g in seq_len(numBio)) {
      
      random_vars <- all.vars(x$random[[g]])
      
      if (length(random_vars) > 1) {
        random_vars <- random_vars[-length(random_vars)]
        biomarker_names <- c(
          paste0("Intercept", g),
          paste0(random_vars, g)
        )
      } else {
        biomarker_names <- paste0("Intercept", g)
      }
      re_names <- c(re_names, biomarker_names)
    }
    
    if (length(re_names) != pREtotal) {
      stop(
        "The number of constructed random-effect names (",
        length(re_names),
        ") does not match the dimension of x$Sig (",
        pREtotal,
        ")."
      )
    }
    
    # ---------------------------------------------------------
    # Calculate standard deviations and correlations
    # ---------------------------------------------------------
    
    if (any(diag(Sig) < 0)) {
      warning(
        "Negative diagonal elements were found in x$Sig; ",
        "the corresponding standard deviations are set to zero."
      )
    }
    
    re_sd <- sqrt(pmax(diag(Sig), 0))
    
    # Calculate correlations safely, including zero-variance cases
    denom <- outer(re_sd, re_sd)
    re_cor <- Sig / denom
    
    re_cor[!is.finite(re_cor)] <- NA_real_
    diag(re_cor) <- ifelse(re_sd > 0, 1, NA_real_)
    
    # ---------------------------------------------------------
    # Construct JM-style display
    # ---------------------------------------------------------
    
    vc <- matrix(
      "",
      nrow = pREtotal,
      ncol = pREtotal
    )
    
    rownames(vc) <- re_names
    
    corr_names <- if (pREtotal > 1) {
      re_names[seq_len(pREtotal - 1)]
    } else {
      character(0)
    }
    
    # Shorter correlation-column labels
    corr_names <- sub("^Intercept", "Intr", corr_names)
    
    colnames(vc) <- c(
      "StdDev",
      corr_names
    )
    
    # Standard deviations
    vc[, 1] <- formatC(
      re_sd,
      format = "f",
      digits = digits
    )
    
    # Lower-triangular correlations
    if (pREtotal > 1) {
      for (i in 2:pREtotal) {
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
