Getinit <- function(cdata, ydata, long.formula, surv.formula,
                    model, ID, RE, survinitial, REML, random, opt) {
  
  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
  cnames <- colnames(cdata)
  ynames <- colnames(ydata)
  
  ##variable check
  if (prod(long %in% ynames) == 0) {
    Fakename <- which(long %in% ynames == FALSE)
    stop(paste0("The variable ", long[Fakename], " not found"))
  }
  if (prod(survival %in% cnames) == 0) {
    Fakename <- which(survival %in% cnames == FALSE)
    stop(paste0("The variable ", survival[Fakename], " not found"))
  }
  if (!(ID %in% ynames)) {
    stop(paste0("ID column ", ID, " not found in the longitudinal dataset!"))
  }
  if (!(ID %in% cnames)) {
    stop(paste0("ID column ", ID, " not found in the survival dataset!"))
  }
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
  
  ## extract covariates
  X <- ydata[, long[-1]]
  X <- as.matrix(cbind(1, X))
  
  ##random effect covariates
  if (model == "interslope") {
    if (prod(RE %in% ynames) == 0) {
      Fakename <- which(RE %in% ynames == FALSE)
      stop(paste0("The variable ", RE[Fakename], " not found"))
    } else {
      #random.xnam <- paste(RE[1:length(RE)], sep = "")
      #fmla.random <- paste("(", paste(random.xnam, collapse= "+"), "|", ID, ")", sep = "")
      p1a <- 1 + length(RE)
      Z <- ydata[, RE]
      Z <- cbind(1, Z)
      Z <- as.matrix(Z)
    }
  } else if (model == "intercept") {
    if (!is.null(RE)) {
      stop("You are fitting a mixed effects location scale model with random intercept only
           but random effects covariates are specified at the same time. Please respecify your model!")
    }
    p1a <- 1
    Z <- rep(1, ydim[1])
    Z <- as.data.frame(Z)
    Z <- as.matrix(Z)
    #fmla.random <- paste("(", 1, "|", ID, ")", sep = "")
    
  } else {
    stop("model should be one of the following options: interslope or intercept.")
  }
  
  Y <- as.vector(ydata[, long[1]])
  X2 <- as.matrix(cdata[, survival[3:length(survival)]])
  survtime <- as.vector(cdata[, survival[1]])
  cmprsk <- as.vector(cdata[, survival[2]])
  
  ## get initial estimates of fixed effects in linear mixed effects model
  # xnam <- paste(long[2:length(long)], sep = "")
  # fmla.fixed <- paste(long[1], " ~ ", paste(xnam, collapse= "+"))
  # fmla <- as.formula(paste(fmla.fixed, fmla.random, sep = "+"))
  # 
  # longfit <- lme4::lmer(fmla, data = ydata, REML = TRUE)
  # beta <- lme4::fixef(longfit)
  # D <- lme4::VarCorr(longfit)
  # name <- names(D)
  # D <- as.data.frame(D[name])
  # D <- as.matrix(D)
  # ##Residuals
  # sigma <- sigma(longfit)^2
  
  if (REML) method <- "REML"
  if (!REML) method <- "ML"
  longfit <- nlme::lme(fixed = long.formula, random = random, data = ydata, method = method, 
                       control = nlme::lmeControl(opt = opt))
  beta <- longfit$coefficients$fixed
  sigma <- longfit$sigma^2
  D <- as.matrix(nlme::getVarCov(longfit))
  
  if (sum(unique(cmprsk)) <= 3) {
    if (prod(c(0, 1, 2) %in% unique(cmprsk))) {
      if (survinitial) {
        ## get initial estimates of fixed effects in competing risks model
        surv_xnam <- paste(survival[3:length(survival)], sep = "")
        survfmla.fixed <- paste(surv_xnam, collapse= "+")
        survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==1)")
        survfmla <- as.formula(paste(survfmla.out1, survfmla.fixed, sep = "~"))
        fitSURV1 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
        gamma1 <- as.vector(fitSURV1$coefficients)
        
        survfmla.out2 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==2)")
        survfmla <- as.formula(paste(survfmla.out2, survfmla.fixed, sep = "~"))
        fitSURV2 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
        gamma2 <- as.vector(fitSURV2$coefficients)  
      } else {
        gamma1 = as.vector(rep(0, ncol(X2)))
        gamma2 = as.vector(rep(0, ncol(X2)))
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
      
      a <- list(beta, gamma1, gamma2, alpha1, alpha2, Sig, sigma, 
                Z, X, Y, X2, survtime, cmprsk, ydata, cdata, mdata)
      
      names(a) <- c("beta", "gamma1", "gamma2", "alpha1", "alpha2", "Sig", "sigma",
                    "Z", "X1", "Y", "X2", "survtime", "cmprsk", "ydata",
                    "cdata", "mdata")
      
      return(a)
      
    } 
    
    if (prod(c(0, 1) %in% unique(cmprsk))) {
      
      if (survinitial) {
        ## get initial estimates of fixed effects in survival model
        surv_xnam <- paste(survival[3:length(survival)], sep = "")
        survfmla.fixed <- paste(surv_xnam, collapse= "+")
        survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], ")")
        survfmla <- as.formula(paste(survfmla.out1, survfmla.fixed, sep = "~"))
        fitSURV1 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
        gamma1 <- as.vector(fitSURV1$coefficients)
      } else {
        gamma1 = as.vector(rep(0, ncol(X2)))
      }
      
      
      if (model == "intercept") {
        alpha1 = as.vector(0)
        Sig <- D
      } else {
        alpha1 = as.vector(rep(0, p1a))
        Sig <- D
      }
      
      a <- list(beta, gamma1, alpha1, Sig, sigma,
                Z, X, Y, X2, survtime, cmprsk, ydata, cdata, mdata)
      
      names(a) <- c("beta", "gamma1", "alpha1", "Sig", "sigma",
                    "Z", "X1", "Y", "X2", "survtime", "cmprsk", "ydata",
                    "cdata", "mdata")
      
      return(a)
      
    }
  } else {
    stop(paste0("The ", survival[2], " variable is specified incorrectly! Program stops."))
  }
  
}