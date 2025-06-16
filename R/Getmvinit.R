Getmvinit <- function(cdata, ydata, long.formula, surv.formula,
                      model, ID, RE, REML, random, opt, initial.para,
                      latAsso = "sre") {
  
  cnames <- colnames(cdata)
  ynames <- colnames(ydata)
  
  survival <- all.vars(surv.formula)
  
  ydim = dim(ydata)
  cdim = dim(cdata)
  
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
  
  for(g in 1:numBio){
    
    ydata <- ydataAll[[g]]
    long <- all.vars(long.formula[[g]])
    
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
      Z <- rep(1, ydim[1])
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
  
  cmprsk <- as.vector(cdata[, survival[2]])
  
  dimmi <- c()
  
  if (sum(unique(cmprsk)) <= 3) {
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
            
            #if(colnames(mi))
            colnames(mi)[1] <- "Intercept"
            colnames(mi) <- paste("random", colnames(mi), g, sep = "_")
            cdata <- cbind(cdata,mi)
            
            mifmla <- paste0(names(cdata)[(length(all.vars(surv.formula))+1 + index):(length(all.vars(surv.formula)) + index + dimmi[g])], collapse = "+")
            index = index + dimmi[g]
            survfmla.comb <- paste0(survfmla.comb, "+", mifmla)
          }
          
          survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==1)")
          survfmla <- as.formula(paste(survfmla.out1, survfmla.comb, sep = "~"))
          fitSURV1 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
          survfmla.out2 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==2)")
          survfmla <- as.formula(paste(survfmla.out2, survfmla.comb, sep = "~"))
          fitSURV2 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
        }
        
        #survfmla <- survival::Surv(survtime, cmprsk == 1) ~X21 + X22 + random_Intercept +
        #random_time
        
        allalphaInd1 <- allalphaInd2 <- c()
        alpha <- list(alpha1 = list(), alpha2 = list())
        
        fitSURV1co <- fitSURV1$coefficients
        fitSURV2co <- fitSURV2$coefficients
        index <- 0
        for (g in 1:numBio) {
          alphaInd1 <- (length(all.vars(surv.formula[[3]]))+index+1):(length(all.vars(surv.formula[[3]]))+index+dimmi[g])
          alphaInd2 <- (length(all.vars(surv.formula[[3]]))+index+1):(length(all.vars(surv.formula[[3]]))+index+dimmi[g])
          alpha[[1]][[g]] <- fitSURV1co[alphaInd1]
          alpha[[2]][[g]] <- fitSURV2co[alphaInd2]
          allalphaInd1 <- c(allalphaInd1, alphaInd1)
          allalphaInd2 <- c(allalphaInd2, alphaInd2)
          index <- index + dimmi[g]
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
            
            #if(colnames(mi))
            colnames(mi)[1] <- "Intercept"
            colnames(mi) <- paste("random", colnames(mi), g, sep = "_")
            cdata <- cbind(cdata,mi)
            
            mifmla <- paste0(names(cdata)[(length(all.vars(surv.formula))+1 + index):(length(all.vars(surv.formula)) + index + dimmi[g])], collapse = "+")
            index = index + dimmi[g]
            survfmla.comb <- paste0(survfmla.comb, "+", mifmla)
          }
          survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==1)")
          survfmla <- as.formula(paste(survfmla.out1, survfmla.comb, sep = "~"))
          fitSURV1 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
          
          allalphaInd1 <- c()
          alpha <- list(alpha1 = list())
          fitSURV1co <- fitSURV1$coefficients
          index <- 0
          for (g in 1:numBio) {
            alphaInd1 <- (length(all.vars(surv.formula[[3]]))+index+1):(length(all.vars(surv.formula[[3]]))+index+dimmi[g])
            alpha[[1]][[g]] <- fitSURV1co[alphaInd1]
            allalphaInd1 <- c(allalphaInd1, alphaInd1)
            index <- index + dimmi[g]
          }
          gamma1 <- fitSURV1co[-allalphaInd1]
          
        }
        
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
                long.formula, surv.formula)
      
      names(a) <- c("beta", "gamma1", "gamma2", "alpha", "Sig", "sigma",
                    "Z", "X1", "Y", "W", "survtime", "cmprsk", "b", "ydata",
                    "cdata", "mdata", "long.formula", "surv.formula")
      
      return(a)
    } else {
      a <- list(beta, gamma1, alpha, Sig, sigma,
                Zlist, X, Y, X2, survtime, cmprsk, bi, ydatanew, cdata, mdata,
                long.formula, surv.formula)
      
      names(a) <- c("beta", "gamma1", "alpha", "Sig", "sigma",
                    "Z", "X1", "Y", "W", "survtime", "cmprsk", "b" , "ydata",
                    "cdata", "mdata", "long.formula", "surv.formula")
      
      return(a)
    }
    
  } else {
    stop(paste0("The current version can only support up to two competing risks events. Program stops."))
  }
  
}