Getinit <- function(cdata, ydata, long.formula, surv.formula,
                    model, ID, RE, survinitial, REML, random, opt) {
  
  cnames <- colnames(cdata)
  ynames <- colnames(ydata)
  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
  
  if (is.null(RE) & model == "interslope") {
    stop("Random effects covariates must be specified.")
  }
  yID <- unique(ydata[, ID])
  cID <- cdata[, ID]
  if (prod(yID == cID) == 0) {
    stop("The order of subjects in ydata doesn't match with cdata.")
  }
  
  ydim = dim(ydata)
  cdim = dim(cdata)
  
  orderdata <- sortdata(cdata, ydata, ID, surv.formula, long.formula)
  
  ydata <- orderdata$ydata
  cdata <- orderdata$cdata
  mdata <- orderdata$mdata
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
      Z <- as.matrix(Z)
    }
  } else if (model == "intercept") {
    if (!is.null(RE)) {
      stop("You are fitting a mixed effects model with random intercept only
           but random effects covariates are specified at the same time. Please respecify your model!")
    }
    p1a <- 1
    Z <- rep(1, ydim[1])
    Z <- as.data.frame(Z)
    Z <- as.matrix(Z)
    
  } else {
    stop("model should be one of the following options: interslope or intercept.")
  }
  
  if (REML) method <- "REML"
  if (!REML) method <- "ML"
  longfit <- nlme::lme(fixed = long.formula, random = random, data = ydata, method = method, 
                       control = nlme::lmeControl(opt = opt))
  beta <- longfit$coefficients$fixed
  sigma <- longfit$sigma^2
  D <- as.matrix(nlme::getVarCov(longfit))
  
  cmprsk <- as.vector(cdata[, survival[2]])
  
  if (sum(unique(cmprsk)) <= 3) {
    if (prod(c(0, 1, 2) %in% unique(cmprsk))) {
      survfmla.fixed <- surv.formula[3]
      survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==1)")
      survfmla <- as.formula(paste(survfmla.out1, survfmla.fixed, sep = "~"))
      fitSURV1 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
      
      survfmla.out2 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==2)")
      survfmla <- as.formula(paste(survfmla.out2, survfmla.fixed, sep = "~"))
      fitSURV2 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
      if (survinitial) {
        gamma1 <- fitSURV1$coefficients
        gamma2 <- fitSURV2$coefficients 
      } else {
        gamma1 = rep(0, length(fitSURV1$coefficients))
        names(gamma1) <- names(fitSURV1$coefficients)
        gamma2 = rep(0, length(fitSURV2$coefficients))
        names(gamma2) <- names(fitSURV2$coefficients)
      }
      
      if (model == "intercept") {
        alpha1 = as.vector(0)
        alpha2 = as.vector(0)
        Sig <- D
        
      } else {
        alpha1 = as.vector(rep(0, p1a))
        alpha2 = as.vector(rep(0, p1a))
        Sig <- D
      }
    } else {
      
      survfmla.fixed <- surv.formula[3]
      survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==1)")
      
      survfmla <- as.formula(paste(survfmla.out1, survfmla.fixed, sep = "~"))
      fitSURV1 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
      if (survinitial) {
        gamma1 <- fitSURV1$coefficients
      } else {
        gamma1 = rep(0, ncol(X2))
        names(gamma1) <- names(fitSURV1$coefficients)
      }
      
      
      if (model == "intercept") {
        alpha1 = as.vector(0)
        Sig <- D
      } else {
        alpha1 = as.vector(rep(0, p1a))
        Sig <- D
      }
      
    } 
    
    getdum <- getdummy(long.formula = long.formula, surv.formula = surv.formula, 
                       random = random, ydata = ydata, cdata = cdata)
    
    ydata <- getdum$ydata
    cdata <- getdum$cdata
    
    ## extract covariates
    X <- ydata[, -c(1:2)]
    X <- as.matrix(cbind(1, X))
    
    Y <- as.vector(ydata[, 2])
    X2 <- as.matrix(cdata[, -c(1:3)])
    survtime <- as.vector(cdata[, survival[1]])
    cmprsk <- as.vector(cdata[, survival[2]])
    
    if (prod(c(0, 1, 2) %in% unique(cmprsk))) {
      a <- list(beta, gamma1, gamma2, alpha1, alpha2, Sig, sigma, 
                Z, X, Y, X2, survtime, cmprsk, ydata, cdata, mdata, 
                long.formula, surv.formula)
      
      names(a) <- c("beta", "gamma1", "gamma2", "alpha1", "alpha2", "Sig", "sigma",
                    "Z", "X1", "Y", "X2", "survtime", "cmprsk", "ydata",
                    "cdata", "mdata", "long.formula", "surv.formula")
      
      return(a)
    } else {
      a <- list(beta, gamma1, alpha1, Sig, sigma,
                Z, X, Y, X2, survtime, cmprsk, ydata, cdata, mdata,
                long.formula, surv.formula)
      
      names(a) <- c("beta", "gamma1", "alpha1", "Sig", "sigma",
                    "Z", "X1", "Y", "X2", "survtime", "cmprsk", "ydata",
                    "cdata", "mdata", "long.formula", "surv.formula")
      
      return(a)
    }
    
  } else {
    stop(paste0("The current version can only support up to two competing risks events. Program stops."))
  }
  
}