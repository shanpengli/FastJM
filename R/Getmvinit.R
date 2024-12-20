Getmvinit <- function(cdata, ydata, long.formula, surv.formula,
                      model, ID, RE, survinitial, REML, random, opt, latAsso = "sre") {
  
  cnames <- colnames(cdata)
  ynames <- colnames(ydata)
  
  survival <- all.vars(surv.formula)
  
  ydim = dim(ydata)
  cdim = dim(cdata)
  
  numBio <- length(long.formula)
  
  beta <- sigma <- D <- bi <- Sig <- mdata <- ydatanew <- list()
  
  
  orderdata <- sortmvdata(cdata, ydata, ID, surv.formula, long.formula)
  ydata <- orderdata$ydata
  cdata <- orderdata$cdata
  # need to change this
  for(i in 1:numBio){
    mdata[[i]] <- orderdata$mdata
  }
  
  Zlist <- list()
  
  for(i in 1:numBio){
    long <- all.vars(long.formula[[i]])
    
    ##random effect covariates
    if (model == "interslope") {
      if (prod(RE %in% ynames) == 0) {
        Fakename <- which(RE %in% ynames == FALSE)
        stop(paste0("The variable ", RE[Fakename], " not found in the longitudinal dataset.\n"))
      } else if (prod(RE %in% long) == 0) {
        Fakename <- which(RE %in% long == FALSE)
        stop(paste0("The variable ", RE[Fakename], " not found in the long.formula argument.
                  Please include this variable in the random argument.\n"))
      } else {
        p1a <- 1 + length(RE)
        Z <- ydata[, RE]
        Z <- cbind(1, Z)
        Zlist[[i]] <- as.matrix(Z)
      }
    } else if (model == "intercept") {
      if (!is.null(RE)) {
        stop("You are fitting a mixed effects model with random intercept only
           but random effects covariates are specified at the same time. Please respecify your model!")
      }
      p1a <- 1
      Z <- rep(1, ydim[1])
      Z <- as.data.frame(Z)
      Zlist[[i]]<- as.matrix(Z)
      
    } else {
      stop("model should be one of the following options: interslope or intercept.")
    }
    
    if (REML) method <- "REML"
    if (!REML) method <- "ML"
    
    
    
    longfit <- try(nlme::lme(fixed = long.formula[[i]], random = random, data = ydata, method = method,
                             control = nlme::lmeControl(opt = opt), na.action = na.omit), silent = TRUE)
    
    if ('try-error' %in% class(longfit)) {
      return(NULL)
    } else {
      beta[[i]] <- longfit$coefficients$fixed
      sigma[[i]] <- longfit$sigma^2
      Sig[[i]] <- as.matrix(nlme::getVarCov(longfit))
      bi[[i]] <- longfit$coefficients$random[[1]]
    }
  
    
    getdum <- getmvdummy(long.formula = long.formula[[i]], surv.formula = surv.formula,
                       random = random, ydata = ydata, cdata = cdata)
    
    ydatanew[[i]] <- getdum$ydata
    cdatanew <- getdum$cdata
    
  }
  
  
  
  cmprsk <- as.vector(cdata[, survival[2]])
  
  dimmi <- c()
  
  if (sum(unique(cmprsk)) <= 3) {
    if (prod(c(0, 1, 2) %in% unique(cmprsk))) {
      
      if(latAsso == "sre"){
        for(i in 1:numBio){
          mi = bi[[i]]
          dimmi[i] <- ncol(mi)
          #if(colnames(mi))
          colnames(mi)[1] <- "Intercept"
          colnames(mi) <- paste("random", colnames(mi), i, sep = "_")
          cdata <- cbind(cdata,mi)
        }
      }
      
      
      
      survfmla.fixed <- surv.formula[3]
      survfmla.fixed <- gsub("\\(+$", "",survfmla.fixed)
      
      
      mifmla <- paste0(names(cdata)[(ncol(cdata)-ncol(mi)*numBio+1):ncol(cdata)], collapse = "+")
      survfmla.fixed <- paste0(survfmla.fixed, "+", mifmla)
      survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==1)")
      survfmla <- as.formula(paste(survfmla.out1, survfmla.fixed, sep = "~"))
      #survfmla <- survival::Surv(survtime, cmprsk == 1) ~X21 + X22 + random_Intercept + 
      #random_time
      fitSURV1 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
      survfmla.out2 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==2)")
      survfmla <- as.formula(paste(survfmla.out2, survfmla.fixed, sep = "~"))
      fitSURV2 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
      allalphaInd1 <- allalphaInd2 <- c()
      alpha <- list(alpha1 = list(), alpha2 = list())
      if (survinitial) {
        fitSURV1co <- fitSURV1$coefficients
        fitSURV2co <- fitSURV2$coefficients
        
        
        for (i in 1:numBio) {
          alphaInd1 <- (length(fitSURV1co)-sum(dimmi[1:i])+1):(length(fitSURV1co)-sum(dimmi[1:i])+dimmi[i])
          alphaInd2 <- (length(fitSURV1co)-sum(dimmi[1:i])+1):(length(fitSURV1co)-sum(dimmi[1:i])+dimmi[i])
          alpha[[numBio-i+1]][[1]] <- fitSURV1co[alphaInd1]
          alpha[[numBio-i+1]][[2]] <- fitSURV2co[alphaInd2]
          allalphaInd1 <- c(allalphaInd1, alphaInd1)
          allalphaInd2 <- c(allalphaInd2, alphaInd2)
        }
        
        gamma1 <- fitSURV1co[-allalphaInd1]
        gamma2 <- fitSURV2co[-allalphaInd2]
      } else {
        gamma1 = rep(0, length(fitSURV1co[-allalphaInd1]))
        names(gamma1) <- names(fitSURV1co)
        gamma2 = rep(0, length(fitSURV2co[-allalphaInd2]))
        names(gamma2) <- names(fitSURV2co)
        
        for(i in 1:numBio){
          alpha[[i]][[1]] <- rep(0, dimmi[i])
          alpha[[i]][[2]] <- rep(0, dimmi[i])
        }
        
      }
      
      
      
      
    } else {
      
      survfmla.fixed <- surv.formula[3]
      survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==1)")
      
      survfmla <- as.formula(paste(survfmla.out1, survfmla.fixed, sep = "~"))
      fitSURV1 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
      if (survinitial) {
        gamma1 <- fitSURV1$coefficients
      } else {
        gamma1 = rep(0, length(fitSURV1$coefficients))
        names(gamma1) <- names(fitSURV1$coefficients)
      }
      alpha <- list()
      if (model == "intercept") {
        alpha1 = as.vector(0)
        # Sig <- D
      } else {
        alpha1 = as.vector(rep(0, p1a))
        # Sig <- D
      }
      
      alpha[[i]] <- list(alpha1)
      
    }
    
    
    
    ## extract covariates
    
    X <- Y <- list()
    
    ## will need to fix this later
    for(i in 1:numBio){
      Xtemp <- ydatanew[[i]][, -c(1,2)] # reorder later
      X[[i]] <- as.matrix(cbind(1, Xtemp))
      Y[[i]] <- as.vector(ydatanew[[i]][, 2])
    }
    
    
    X2 <- as.matrix(cdatanew[, -c(1:3)])
    survtime <- as.vector(cdatanew[, survival[1]])
    cmprsk <- as.vector(cdatanew[, survival[2]])
    
    if (prod(c(0, 1, 2) %in% unique(cmprsk))) {
      a <- list(beta, gamma1, gamma2, alpha, Sig, sigma,
                Zlist, X, Y, X2, survtime, cmprsk, bi, ydata, cdata, mdata,
                long.formula, surv.formula)
      
      names(a) <- c("beta", "gamma1", "gamma2", "alpha", "Sig", "sigma",
                    "Z", "X1", "Y", "W", "survtime", "cmprsk", "b", "ydata",
                    "cdata", "mdata", "long.formula", "surv.formula")
      
      return(a)
    } else {
      a <- list(beta, gamma1, alpha1, Sig, sigma,
                Zlist, X, Y, X2, survtime, cmprsk, bi, ydata, cdata, mdata,
                long.formula, surv.formula)
      
      names(a) <- c("beta", "gamma1", "alpha1", "Sig", "sigma",
                    "Z", "X1", "Y", "X2", "survtime", "cmprsk", "b" , "ydata",
                    "cdata", "mdata", "long.formula", "surv.formula")
      
      return(a)
    }
    
  } else {
    stop(paste0("The current version can only support up to two competing risks events. Program stops."))
  }
  
}