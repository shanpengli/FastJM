Getmvinit <- function(cdata, ydata, long.formula, surv.formula,
                        model, ID, RE, REML, random, opt, initial.para,
                        latAsso = "sre", landmark = TRUE,  s = NULL, ytime = NULL) {
  
  survival <- all.vars(surv.formula)
  
  numBio <- length(long.formula)
  
  mdata <- ydatanew <- bi <- list()
  
  if (is.null(initial.para)) {
    beta <- D <- Sig <- list()
    sigma <- c()
  } else {
    beta <- initial.para$beta
    sigma <- initial.para$sigma
    Sig <- initial.para$Sig
  }
  
  orderdata <- sortmvdata(cdata, ydata, ID, surv.formula, long.formula)
  ydataAll <- orderdata$ydata
  cdata <- orderdata$cdata
  mdata <- orderdata$mdata
  
  Zlist <- list()
  
  if (latAsso %in% c("present", "presentlp") && !landmark) {
    stop("latAsso = 'present' or 'presentlp' requires landmark = TRUE and a valid s.")
  }
  
  if (latAsso == "present" && is.null(s)) {
    stop("latAsso = 'present' requires s so X(s) can be constructed.")
  }
  
  
  if (landmark) {
    
    if (is.null(s)) {
      warning("landmark = TRUE but no landmark time provided. Landmark preprocessing is skipped.")
    }
    else{
      
      if (is.null(ytime)) {
        stop("For landmarking, you must provide the longitudinal time variable name through 'ytime'.")
      }
      
      survtime.var <- survival[1]
      status.var <- survival[2]
      
      # get subjects in after landmark
      
      keep.id <- cdata[[ID]][cdata[[survival[1]]] > s]
      keep.id <- unique(keep.id)
      
      cdata <- cdata[cdata[[ID]] %in% keep.id, , drop = FALSE]
      
      # restrict data to those who are at risk
      ydataAll <- lapply(ydataAll, function(dat) {
        dat <- dat[dat[[ID]] %in% keep.id, , drop = FALSE]
        dat
      })
      
      # # get ids unique
      ids.with.long <- unique(unlist(lapply(ydataAll, function(dat) unique(dat[[ID]]))))
      cdata <- cdata[cdata[[ID]] %in% ids.with.long, , drop = FALSE]
      
    }
    
  }
  
  
  
  cdim = dim(cdata)
  cnames <- colnames(cdata)
  
  
  for(g in 1:numBio){
    
    ydata <- ydataAll[[g]]
    long <- all.vars(long.formula[[g]])
    
    #  = dim(ydata)
    ynames <- colnames(ydata)
    
    
    
    
    ##random effect covariates
    if (model[[g]] == "interslope") {
      if (prod(RE[[g]] %in% ynames) == 0) {
        Fakename <- which(RE[[g]] %in% ynames == FALSE)
        stop(paste0("The variable ", RE[[g]][Fakename], " not found in the longitudinal dataset.\n"))
      } else if (prod(RE[[g]] %in% long) == 0) {
        Fakename <- which(RE[[g]] %in% long == FALSE)
        stop(paste0("The variable ", RE[[g]][Fakename], " not found in the long.formula argument.
                  Please include this variable in the random argument.\n"))
      } else {
        p1a <- 1 + length(RE[[g]])
        Z <- ydata[, RE[[g]]]
        Z <- cbind(1, Z)
        Zlist[[g]] <- as.matrix(Z)
      }
    } else if (model[[g]] == "intercept") {
      if (!is.null(RE[[g]])) {
        stop("You are fitting a mixed effects model with random intercept only
           but random effects covariates are specified at the same time. Please respecify your model!")
      }
      p1a <- 1
      Z <- rep(1, nrow(ydata))
      Z <- as.data.frame(Z)
      Zlist[[g]]<- as.matrix(Z)
      
    } else {
      stop("model should be one of the following options: interslope or intercept.")
    }
    
    if (REML) method <- "REML"
    if (!REML) method <- "ML"
    
    if (is.null(initial.para)) {
      
      longfit <- try(nlme::lme(fixed = long.formula[[g]], random = random[[g]], data = ydata, method = method,
                               control = nlme::lmeControl(opt = opt), na.action = na.omit), silent = TRUE)
      
      if ('try-error' %in% class(longfit)) {
        return(NULL)
      } else {
        beta[[g]] <- longfit$coefficients$fixed
        sigma[g] <- longfit$sigma^2
        Sig[[g]] <- as.matrix(nlme::getVarCov(longfit))
        bi[[g]] <- longfit$coefficients$random[[1]]
      }
      
    }
    
    getdum <- getmvdummy(long.formula = long.formula[[g]], surv.formula = surv.formula,
                         random = random[[g]], ydata = ydata, cdata = cdata)
    
    ydatanew[[g]] <- getdum$ydata
    cdatanew <- getdum$cdata
    
    
  }
  
  mdata <- lapply(ydatanew, function(dat) {
    counts <- table(dat[[ID]])
    data.frame(
      ID = names(counts),
      ni = as.integer(counts),
      row.names = NULL
    )
  })
  
  for (g in seq_along(mdata)) {
    idx <- match(cdata[[ID]], mdata[[g]][[ID]])
    if (any(is.na(idx))) {
      stop(paste("Subjects in cdata missing from ydatanew for biomarker", g))
    }
    mdata[[g]] <- mdata[[g]][idx, , drop = FALSE]
  }
  
  cmprsk <- as.vector(cdata[, survival[2]])
  
  dimmi <- c()
  
  Zs<- NULL
  
  if (landmark && latAsso %in% c("present", "presentlp")) {
    
    pREvec <- sapply(Zlist, ncol)
    pREtotal <- sum(pREvec)
    
    Zs <- matrix(0, nrow = numBio, ncol = pREtotal)
    
    index <- 0
    for (g in seq_len(numBio)) {
      qg <- pREvec[g]
      
      Zg_s <- rep(0, qg)
      
      if (qg >= 1) Zg_s[1] <- 1
      if (qg >= 2) Zg_s[2] <- s
      if (qg >= 3) Zg_s[3] <- s^2
      
      Zs[g, (index + 1):(index + qg)] <- Zg_s
      
      index <- index + qg
    }
  }
  Xs <- NULL
  
  if (landmark && latAsso == "present") {
    
    Xs <- vector("list", numBio)
    
    for (g in seq_len(numBio)) {
      
      dat_g <- ydatanew[[g]]
      dat_g <- dat_g[order(dat_g[[ID]], dat_g[[ytime]]), , drop = FALSE]
      
      last_rows <- dat_g[!duplicated(dat_g[[ID]], fromLast = TRUE), , drop = FALSE]
      last_rows <- last_rows[match(cdata[[ID]], last_rows[[ID]]), , drop = FALSE]
      
      if (any(is.na(last_rows[[ID]]))) {
        stop("Some subjects in cdata are missing longitudinal data for biomarker ", g)
      }
      
      last_rows[[ytime]] <- s
      
      Xs[[g]] <- model.matrix(
        delete.response(terms(long.formula[[g]])),
        data = last_rows
      )
    }
  }
  
  if (all(unique(cmprsk) %in% c(0, 1, 2))) {
    if (prod(c(0, 1, 2) %in% unique(cmprsk))) {
      
      if (is.null(initial.para)) {
        
        survfmla.fixed <- surv.formula[3]
        survfmla.fixed <- gsub("\\(+$", "",survfmla.fixed)
        survfmla.comb <- survfmla.fixed
        
        if(latAsso == "sre"){
          
          index = 1
          for(g in 1:numBio){
            mi = bi[[g]]
            dimmi[g] <- ncol(mi)
            
            colnames(mi)[1] <- "Intercept"
            colnames(mi) <- paste("random", colnames(mi), g, sep = "_")
            cdata <- cbind(cdata, mi)
            
            mifmla <- paste0(
              names(cdata)[(length(all.vars(surv.formula)) + 1 + index):
                             (length(all.vars(surv.formula)) + index + dimmi[g])],
              collapse = "+"
            )
            
            index = index + dimmi[g]
            survfmla.comb <- paste0(survfmla.comb, "+", mifmla)
          }
          
        } else if(latAsso %in% c("present", "presentlp")){
          
          Binit <- do.call(cbind, lapply(seq_len(nrow(cdata)), function(i) {
            unlist(lapply(seq_len(numBio), function(g) {
              as.numeric(bi[[g]][i, ])
            }))
          }))
          
          randomMat <- Zs %*% Binit
          
          assocInit <- matrix(0, nrow = nrow(cdata), ncol = numBio)
          
          if(latAsso == "presentlp"){
            
            assocInit <- t(randomMat)
            
          } else if(latAsso == "present"){
            
            betaFull <- unlist(beta)
            
            for(i in seq_len(nrow(cdata))){
              fixed_i <- numeric(numBio)
              
              betaIndex <- 0
              for(g in seq_len(numBio)){
                pg <- ncol(Xs[[g]])
                beta_g <- betaFull[(betaIndex + 1):(betaIndex + pg)]
                
                fixed_i[g] <- sum(Xs[[g]][i, ] * beta_g)
                
                betaIndex <- betaIndex + pg
              }
              
              assocInit[i, ] <- fixed_i + randomMat[, i]
            }
          }
          
          index = 1
          for(g in seq_len(numBio)){
            dimmi[g] <- 1
            
            cname <- paste(latAsso, "bio", g, sep = "_")
            cdata[[cname]] <- assocInit[, g]
            
            mifmla <- cname
            index = index + dimmi[g]
            
            survfmla.comb <- paste0(survfmla.comb, "+", mifmla)
          }
        }
        
        survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==1)")
        survfmla <- as.formula(paste(survfmla.out1, survfmla.comb, sep = "~"))
        fitSURV1 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
        
        survfmla.out2 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==2)")
        survfmla <- as.formula(paste(survfmla.out2, survfmla.comb, sep = "~"))
        fitSURV2 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
        
        #survfmla <- survival::Surv(survtime, cmprsk == 1) ~X21 + X22 + random_Intercept +
        #random_time
        
        # for (g in 1:numBio) {
        #   alphaInd1 <- (length(all.vars(surv.formula[[3]]))+index+1):(length(all.vars(surv.formula[[3]]))+index+dimmi[g])
        #   alphaInd2 <- (length(all.vars(surv.formula[[3]]))+index+1):(length(all.vars(surv.formula[[3]]))+index+dimmi[g])
        #   alpha[[1]][[g]] <- fitSURV1co[alphaInd1]
        #   alpha[[2]][[g]] <- fitSURV2co[alphaInd2]
        #   allalphaInd1 <- c(allalphaInd1, alphaInd1)
        #   allalphaInd2 <- c(allalphaInd2, alphaInd2)
        #   index <- index + dimmi[g]
        # }
        
        
        fitSURV1co <- fitSURV1$coefficients
        fitSURV2co <- fitSURV2$coefficients
        
        allalphaInd1 <- allalphaInd2 <- c()
        alpha <- list(alpha1 = list(), alpha2 = list())
        
        index <- 0
        base_p <- length(all.vars(surv.formula[[3]]))
        
        for (g in seq_len(numBio)) {
          
          if (latAsso == "sre") {
            qg <- dimmi[g]       # vector alpha per biomarker
          } else {
            qg <- 1              # scalar alpha per biomarker
          }
          
          alphaInd1 <- (base_p + index + 1):(base_p + index + qg)
          alphaInd2 <- (base_p + index + 1):(base_p + index + qg)
          
          alpha[[1]][[g]] <- fitSURV1co[alphaInd1]
          alpha[[2]][[g]] <- fitSURV2co[alphaInd2]
          
          allalphaInd1 <- c(allalphaInd1, alphaInd1)
          allalphaInd2 <- c(allalphaInd2, alphaInd2)
          
          index <- index + qg
        }
        
        gamma1 <- fitSURV1co[-allalphaInd1]
        gamma2 <- fitSURV2co[-allalphaInd2]
        
      } else {
        gamma1 <- initial.para$gamma1
        gamma2 <- initial.para$gamma2
        alpha <- initial.para$alpha
        
      }
      
    } else {
      
      
      if (is.null(initial.para)) {
        
        survfmla.fixed <- surv.formula[3]
        survfmla.fixed <- gsub("\\(+$", "",survfmla.fixed)
        survfmla.comb <- survfmla.fixed
        
        allalphaInd1 <- c()
        alpha <- list(alpha1 = list())
        
        
        if(latAsso == "sre"){
          
          index = 1
          for(g in 1:numBio){
            mi = bi[[g]]
            dimmi[g] <- ncol(mi)
            
            colnames(mi)[1] <- "Intercept"
            colnames(mi) <- paste("random", colnames(mi), g, sep = "_")
            cdata <- cbind(cdata, mi)
            
            mifmla <- paste0(
              names(cdata)[(length(all.vars(surv.formula)) + 1 + index):
                             (length(all.vars(surv.formula)) + index + dimmi[g])],
              collapse = "+"
            )
            
            index = index + dimmi[g]
            survfmla.comb <- paste0(survfmla.comb, "+", mifmla)
          }
          
        } else if(latAsso %in% c("present", "presentlp")){
          
          ## Build initial association value per subject and biomarker
          ## presentlp: Z(s)b
          ## present:   X(s)beta + Z(s)b
          
          Binit <- do.call(cbind, lapply(seq_len(nrow(cdata)), function(i) {
            unlist(lapply(seq_len(numBio), function(g) {
              as.numeric(bi[[g]][i, ])
            }))
          }))
          
          randomMat <- Zs %*% Binit   # numBio x numSubj
          
          assocInit <- matrix(0, nrow = nrow(cdata), ncol = numBio)
          
          if(latAsso == "presentlp"){
            
            assocInit <- t(randomMat)
            
          } else if(latAsso == "present"){
            
            betaFull <- unlist(beta)
            
            for(i in seq_len(nrow(cdata))){
              fixed_i <- numeric(numBio)
              
              betaIndex <- 0
              for(g in seq_len(numBio)){
                pg <- ncol(Xs[[g]])
                beta_g <- betaFull[(betaIndex + 1):(betaIndex + pg)]
                
                fixed_i[g] <- sum(Xs[[g]][i, ] * beta_g)
                
                betaIndex <- betaIndex + pg
              }
              
              assocInit[i, ] <- fixed_i + randomMat[, i]
            }
          }
          
          index = 1
          for(g in 1:numBio){
            dimmi[g] <- 1
            
            cname <- paste(latAsso, "bio", g, sep = "_")
            cdata[[cname]] <- assocInit[, g]
            
            mifmla <- cname
            index = index + dimmi[g]
            
            survfmla.comb <- paste0(survfmla.comb, "+", mifmla)
          }
        }
        
        survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==1)")
        survfmla <- as.formula(paste(survfmla.out1, survfmla.comb, sep = "~"))
        fitSURV1 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
        
        fitSURV1co <- fitSURV1$coefficients
        
        allalphaInd1 <- c()
        alpha <- list(alpha1 = list())
        
        index <- 0
        base_p <- length(all.vars(surv.formula[[3]]))
        
        for (g in seq_len(numBio)) {
          
          if (latAsso == "sre") {
            qg <- dimmi[g]
          } else {
            qg <- 1
          }
          
          alphaInd1 <- (base_p + index + 1):(base_p + index + qg)
          
          alpha[[1]][[g]] <- fitSURV1co[alphaInd1]
          allalphaInd1 <- c(allalphaInd1, alphaInd1)
          
          index <- index + qg
        }
        
        gamma1 <- fitSURV1co[-allalphaInd1]
        
        
      } else {
        gamma1 <- initial.para$gamma1
        alpha <- initial.para$alpha
      }
    }
    
    
    ## extract covariates
    X <- Y <- list()
    
    ## will need to fix this later
    for(g in 1:numBio){
      Xtemp <- ydatanew[[g]][, -c(1,2)] # reorder later
      X[[g]] <- as.matrix(cbind(1, Xtemp))
      Y[[g]] <- as.vector(ydatanew[[g]][, 2])
    }
    
    
    X2 <- as.matrix(cdatanew[, -c(1:3)]) # prob want to check this part again
    survtime <- as.vector(cdatanew[, survival[1]])
    cmprsk <- as.vector(cdatanew[, survival[2]])
    
    
    
    if (prod(c(0, 1, 2) %in% unique(cmprsk))) {
      
      
      a <- list(beta, gamma1, gamma2, alpha, Sig, sigma,
                Zlist, X, Y, X2, survtime, cmprsk, bi, ydatanew, cdata, mdata,
                long.formula, surv.formula, Zs, Xs)
      
      names(a) <- c("beta", "gamma1", "gamma2", "alpha", "Sig", "sigma",
                    "Z", "X1", "Y", "W", "survtime", "cmprsk", "b", "ydata",
                    "cdata", "mdata", "long.formula", "surv.formula", "Zs", "Xs")
      
      return(a)
    } else {
      a <- list(beta, gamma1, alpha, Sig, sigma,
                Zlist, X, Y, X2, survtime, cmprsk, bi, ydatanew, cdata, mdata,
                long.formula, surv.formula, Zs, Xs)
      
      names(a) <- c("beta", "gamma1", "alpha", "Sig", "sigma",
                    "Z", "X1", "Y", "W", "survtime", "cmprsk", "b" , "ydata",
                    "cdata", "mdata", "long.formula", "surv.formula", "Zs", "Xs")
      
      return(a)
    }
    
  } else {
    stop(paste0("The current version can only support up to two competing risks events. Program stops."))
  }
  
}