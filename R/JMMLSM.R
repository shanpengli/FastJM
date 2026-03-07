##' Joint modeling of longitudinal continuous data and competing risks
##' @title Joint Modeling for Continuous outcomes
##' @param ydata a longitudinal data frame in long format.
##' @param cdata a survival data frame with competing risks or single failure.
##' Each subject has one data entry.
##' @param long.formula a formula object with the response variable and fixed effects covariates
##' to be included in the longitudinal sub-model.
##' @param surv.formula a formula object with the survival time, event indicator, and the covariates
##' to be included in the survival sub-model.
##' @param variance.formula an one-sided formula object with the fixed effects covariates to model the variance of longitudinal sub-model.
##' @param random a one-sided formula object describing the random effects part of the longitudinal sub-model.
##' For example, fitting a random intercept model takes the form ~ 1|ID.
##' Alternatively. Fitting a random intercept and slope model takes the form ~ x1 + ... + xn|ID.
##' @param maxiter the maximum number of iterations of the EM algorithm that the function will perform. Default is 10000.
##' @param epsilon Tolerance parameter. Default is 0.0001.
##' @param quadpoint the number of Gauss-Hermite quadrature points
##' to be chosen for numerical integration. Default is 15 which produces stable estimates in most dataframes.
##' @param print.para Print detailed information of each iteration. Default is FALSE, i.e., not to print the iteration details.
##' @param survinitial Fit a Cox model to obtain initial values of the parameter estimates. Default is TRUE.
##' @param initial.para a list of initialized parameters for EM iteration. Default is NULL. 
##' @param method Method for proceeding numerical integration in the E-step. Default is adaptive.
##' @param opt Optimization method to fit a linear mixed effects model, either nlminb (default) or optim.
##' @param initial.optimizer Method for numerical optimization to be used. Default is \code{BFGS}. 
##' @return  Object of class \code{JMMLSM} with elements
##' \item{ydata}{the input longitudinal dataset for fitting a joint model.
##' It has been re-ordered in accordance with descending observation times in \code{cdata}.}
##' \item{cdata}{the input survival dataset for fitting a joint model.
##' It has been re-ordered in accordance with descending observation times.}
##' \item{PropEventType}{a frequency table of number of events.}
##' \item{beta}{the vector of fixed effects for the mean trajectory in the mixed effects location and scale model.} 
##' \item{tau}{the vector of fixed effects for the within-subject variability in the mixed effects location and scale model.} 
##' \item{gamma1}{the vector of fixed effects for type 1 failure for the survival model.}
##' \item{gamma2}{the vector of fixed effects for type 2 failure for the survival model. 
##' Valid only if \code{CompetingRisk = TRUE}.}
##' \item{alpha1}{the vector of association parameter(s) for the mean trajectory for type 1 failure.}
##' \item{alpha2}{the vector of association parameter(s) for the mean trajectory for type 2 failure. Valid only if \code{CompetingRisk = TRUE}.}
##' \item{vee1}{the vector of association parameter(s) for the within-subject variability for type 1 failure.}
##' \item{vee2}{the vector of association parameter(s) for the within-subject variability for type 2 failure. Valid only if \code{CompetingRisk = TRUE}.}
##' \item{H01}{the matrix that collects baseline hazards evaluated at each uncensored event time for type 1 failure. 
##' The first column denotes uncensored event times, the second column the number of events, and the third columns 
##' the hazards obtained by Breslow estimator.}
##' \item{H02}{the matrix that collects baseline hazards evaluated at each uncensored event time for type 2 failure. 
##' The data structure is the same as \code{H01}. Valid only if \code{CompetingRisk = TRUE}.}
##' \item{Sig}{the variance-covariance matrix of the random effects.}
##' \item{iter}{the total number of iterations until convergence.}
##' \item{convergence}{convergence identifier: 1 corresponds to successful convergence, 
##' whereas 0 to a problem (i.e., when 0, usually more iterations are required).}
##' \item{vcov}{the variance-covariance matrix of all the fixed effects for both models.}
##' \item{sebeta}{the standard error of \code{beta}.}
##' \item{setau}{the standard error of \code{tau}.}
##' \item{segamma1}{the standard error of \code{gamma1}.}
##' \item{segamma2}{the standard error of \code{gamma2}. 
##' Valid only if \code{CompetingRisk = TRUE}.}
##' \item{sealpha1}{the standard error of \code{alpha1}.}
##' \item{sealpha2}{the standard error of \code{alpha2}. Valid only if \code{CompetingRisk = TRUE}.}
##' \item{sevee1}{the standard error of \code{vee1}.}
##' \item{sevee2}{the standard error of \code{vee2}. Valid only if \code{CompetingRisk = TRUE}.}
##' \item{seSig}{the vector of standard errors of covariance of random effects.}
##' \item{loglike}{the log-likelihood value.}
##' \item{EFuntheta}{a list with the expected values of all the functions of random effects.}
##' \item{CompetingRisk}{logical value; TRUE if a competing event are accounted for.}
##' \item{quadpoint}{the number of Gauss Hermite quadrature points used for numerical integration.}
##' \item{LongitudinalSubmodelmean}{the component of the \code{long.formula}.}
##' \item{LongitudinalSubmodelvariance}{the component of the \code{variance.formula}.}
##' \item{SurvivalSubmodel}{the component of the \code{surv.formula}.}
##' \item{random}{the component of the \code{random}.}
##' \item{call}{the matched call.}
##' @examples
##' require(FastJM)
##' data(ydata)
##' data(cdata)
##' ## fit a joint model
##' \dontrun{
##' fit <- JMMLSM(cdata = cdata, ydata = ydata, 
##'               long.formula = Y ~ Z1 + Z2 + Z3 + time,
##'               surv.formula = Surv(survtime, cmprsk) ~ var1 + var2 + var3,
##'               variance.formula = ~ Z1 + Z2 + Z3 + time, 
##'               quadpoint = 6, random = ~ 1|ID, print.para = FALSE)
##'               
##' ## make dynamic prediction of two subjects
##' cnewdata <- cdata[cdata$ID %in% c(122, 152), ]
##' ynewdata <- ydata[ydata$ID %in% c(122, 152), ]
##' survfit <- survfitJMMLSM(fit, seed = 100, ynewdata = ynewdata, cnewdata = cnewdata, 
##'                          u = seq(5.2, 7.2, by = 0.5), Last.time = "survtime",
##'                          obs.time = "time", method = "GH")
##' oldpar <- par(mfrow = c(2, 2), mar = c(5, 4, 4, 4))
##' plot(survfit, include.y = TRUE)
##' par(oldpar)
##' }
##' @export
##' 

JMMLSM <- function(cdata, ydata,
                   long.formula,
                   surv.formula,
                   variance.formula, 
                   random,
                   maxiter = 1000, epsilon = 1e-04, 
                   quadpoint = NULL, print.para = FALSE,
                   survinitial = TRUE,
                   initial.para = NULL,
                   method = "adaptive",
                   opt = "nlminb",
                   initial.optimizer = "BFGS") {
  
  
  if (!inherits(long.formula, "formula") || length(long.formula) != 3) {
    stop("\nMean sub-part of location scale model must be a formula of the form \"resp ~ pred\"")
  }
  if (!inherits(variance.formula, "formula") || length(variance.formula) != 2) {
    stop("\nVariance sub-part of location scale model must be a formula of the form \" ~ pred\"")
  }
  
  if (!inherits(surv.formula, "formula") || length(surv.formula) != 3) {
    stop("\nCox proportional hazards model must be a formula of the form \"Surv(.,.) ~ pred\"")
  }
  
  if (method == "adaptive" & is.null(quadpoint)) {
    quadpoint <- 6
  }
  
  if (method == "standard" & is.null(quadpoint)) {
    quadpoint <- 20
  }
  
  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
  variance <- all.vars(variance.formula)
  random.form <- all.vars(random)
  ID <- random.form[length(random.form)]
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
  if (prod(variance %in% ynames) == 0) {
    Fakename <- which(variance %in% ynames == FALSE)
    stop(paste0("The within-subject variables ", long[Fakename], " not found"))
  }
  if (!(ID %in% ynames)) {
    stop(paste0("ID column ", ID, " not found in the longitudinal dataset!"))
  }
  if (!(ID %in% cnames)) {
    stop(paste0("ID column ", ID, " not found in the survival dataset!"))
  }
  
  
  if (length(random.form) == 1) {
    RE <- NULL
    model <- "intercept"
  } else {
    RE <- random.form[-length(random.form)]
    model <- "interslope"
  }
  
  longfmla <- long.formula
  survfmla <- surv.formula
  varformula <- variance.formula
  rawydata <- ydata
  rawcdata <- cdata
  
  getinit <- Getinit.JMH(cdata = cdata, ydata = ydata, long.formula = long.formula,
                     surv.formula = surv.formula, variance.formula = variance.formula,
                     model = model, ID = ID, RE = RE, random = random, survinitial = survinitial,
                     initial.para = initial.para, opt = opt)

  
  cdata <- getinit$cdata
  ydata <- getinit$ydata
  
  survival <- all.vars(surv.formula)
  status <- as.vector(cdata[, survival[2]])
  
  if (prod(c(0, 1, 2) %in% unique(status))) {
    ## initialize parameters
    
    getriskset <- Getriskset(cdata = cdata, surv.formula = surv.formula)
    
    ## number of distinct survival time
    H01 <- getriskset$tablerisk1
    H02 <- getriskset$tablerisk2
    
    ## initialize parameters
    beta <- getinit$beta
    namesbeta <- names(beta)
    tau <- getinit$tau
    namestau <- names(tau)
    gamma1 <- getinit$gamma1
    namesgamma1 <- names(gamma1)
    gamma2 <- getinit$gamma2
    alpha1 <- getinit$alpha1
    alpha2 <- getinit$alpha2
    vee1 <- getinit$vee1
    vee2 <- getinit$vee2
    Sig <- getinit$Sig
    p1a <- ncol(Sig) - 1
    
    CompetingRisk <- TRUE
  } else {
    ## initialize parameters
    
    getriskset <- Getriskset(cdata = cdata, surv.formula = surv.formula)
    
    ## number of distinct survival time
    H01 <- as.matrix(getriskset$tablerisk1)
    
    ## initialize parameters
    beta <- getinit$beta
    namesbeta <- names(beta)
    tau <- getinit$tau
    namestau <- names(tau)
    gamma1 <- getinit$gamma1
    namesgamma1 <- names(gamma1)
    alpha1 <- getinit$alpha1
    vee1 <- getinit$vee1
    Sig <- getinit$Sig
    p1a <- ncol(Sig) - 1
    
    CompetingRisk <- FALSE
  }
  
  getGH <- GetGHmatrix.JMH(quadpoint = quadpoint, p1a = p1a)
  
  xsmatrix <- getGH$xsmatrix
  wsmatrix <- getGH$wsmatrix
  
  iter=0
  
  Z <- getinit$Z
  X1 <- getinit$X1
  W <- getinit$W
  Y <- getinit$Y
  X2 <- getinit$X2
  survtime <- getinit$survtime
  cmprsk <- getinit$cmprsk
  mdata <- getinit$mdata
  n <- nrow(mdata)
  mdata <- as.data.frame(mdata)
  mdata <- as.vector(mdata$ni)
  mdataS <- rep(0, n) 
  mdataS[1] <- 1
  mdataCum <- cumsum(mdata)
  mdata2 <- mdata - 1
  mdataS[2:n] <- mdataCum[2:n] - mdata2[2:n]
  
  if (CompetingRisk == TRUE) {
    repeat
    {
      iter <- iter + 1
      prebeta <- beta
      pretau <- tau
      pregamma1 <- gamma1
      pregamma2 <- gamma2
      prealpha1 <- alpha1
      prealpha2 <- alpha2
      prevee1 <- vee1
      prevee2 <- vee2
      preSig <- Sig
      preH01 <- H01
      preH02 <- H02
      if (print.para) {
        writeLines("iter is:")
        print(iter)
        writeLines("beta is:")
        print(beta)
        writeLines("tau is:")
        print(tau)
        writeLines("gamma1 is:")
        print(gamma1)
        writeLines("gamma2 is:")
        print(gamma2)
        writeLines("alpha1 is:")
        print(alpha1)
        writeLines("alpha2 is:")
        print(alpha2)
        writeLines("vee1 is:")
        print(vee1)
        writeLines("vee2 is:")
        print(vee2)
        writeLines("Sig is:")
        print(Sig)
      }
      
      if (method == "standard") {
        GetEfun <- GetE.JMH(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02,
                        Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix)
      } else if (method == "adaptive") {
        GetEfun <- GetEad.JMH(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02,
                           Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, initial.optimizer)
      } else {
        stop("Please choose one of the following methods for numerical integration in the E-step: standard, adaptive.")
      }

      GetMpara <- GetM.JMH(GetEfun, beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, 
                       Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS)
      
      beta <- GetMpara$beta
      tau <- GetMpara$tau
      gamma1 <- GetMpara$gamma1
      gamma2 <- GetMpara$gamma2
      alpha1 <- GetMpara$alpha1
      alpha2 <- GetMpara$alpha2
      vee1 <- GetMpara$vee1
      vee2 <- GetMpara$vee2
      Sig <- GetMpara$Sig
      H01 <- GetMpara$H01
      H02 <- GetMpara$H02
      
      if((Diffrelative.JMH(beta, prebeta, tau, pretau, gamma1, pregamma1, gamma2, pregamma2,
               alpha1, prealpha1, alpha2, prealpha2, vee1, prevee1, vee2, prevee2,
               Sig, preSig, H01, preH01, H02, preH02, epsilon) == 0) || (iter == maxiter) || (!is.list(GetEfun))
         || (!is.list(GetMpara))) {
        break
      }
      
    }
    
    if (iter == maxiter) {
      writeLines("program stops because of nonconvergence")
      convergence = 0
      result <- list(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, 
                     H02, Sig, iter, convergence)
      
      names(result) <- c("beta", "tau", "gamma1", "gamma2", "alpha1", "alpha2", "vee1",
                         "vee2", "H01", "H02", "Sig", "iter", "convergence")
      
      class(result) <- "JMMLSM"
      
      return(result)
    } else if (!is.list(GetEfun)) {
      writeLines("Something wrong in the E steps")
      beta <- NULL
      tau <- NULL
      gamma1 <- NULL
      gamma2 <- NULL
      alpha1 <- NULL
      alpha2 <- NULL
      vee1 <- NULL
      vee2 <- NULL
      H01 <- NULL
      H02 <- NULL
      Sig <- NULL
      iter <- NULL
      result <- list(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, 
                     H02, Sig, iter)
      names(result) <- c("beta", "tau", "gamma1", "gamma2", "alpha1", "alpha2", "vee1",
                         "vee2", "H01", "H02", "Sig", "iter")
      
      class(result) <- "JMMLSM"
      
      return(result)
    } else if (!is.list(GetMpara)) {
      writeLines("Something wrong in the M steps")
      beta <- NULL
      tau <- NULL
      gamma1 <- NULL
      gamma2 <- NULL
      alpha1 <- NULL
      alpha2 <- NULL
      vee1 <- NULL
      vee2 <- NULL
      H01 <- NULL
      H02 <- NULL
      Sig <- NULL
      iter <- NULL
      result <- list(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, 
                     H02, Sig, iter)
      names(result) <- c("beta", "tau", "gamma1", "gamma2", "alpha1", "alpha2", "vee1",
                         "vee2", "H01", "H02", "Sig", "iter")
      
      class(result) <- "JMMLSM"
      
      return(result)
    } else {
      
      convergence = 1
      
      if (method == "standard") {
        GetEfun <- GetE.JMH(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02,
                        Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix)
      } else if (method == "adaptive") {
        GetEfun <- GetEad.JMH(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02,
                          Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, initial.optimizer)
      } else {
        stop("Please choose one of the following methods for numerical integration in the E-step: standard, adaptive.")
      }
      
      FUNENW <- as.vector(GetEfun$FUNENW)
      FUNBENW <- as.matrix(GetEfun$FUNBENW)
      FUNBS <- as.matrix(GetEfun$FUNBS)
      FUNBW <- as.matrix(GetEfun$FUNBW)
      FUNWS <- as.vector(GetEfun$FUNWS)
      FUNBSENW <- as.matrix(GetEfun$FUNBSENW) 
      FUNEC <- as.matrix(GetEfun$FUNEC)
      FUNBEC <- as.matrix(GetEfun$FUNBEC)
      FUNBSEC <- as.matrix(GetEfun$FUNBSEC)
      FUNWEC <- as.matrix(GetEfun$FUNWEC)
      FUNWSEC <- as.matrix(GetEfun$FUNWSEC)
      FUNB <- as.matrix(GetEfun$FUNB)
      FUNW <- as.vector(GetEfun$FUNW)
      
      EFuntheta <- list(FUNENW = FUNENW, FUNBENW = FUNBENW, FUNBS = FUNBS,
                        FUNBW = FUNBW, FUNWS = FUNWS, FUNBSENW = FUNBSENW,
                        FUNEC = FUNEC, FUNBEC = FUNBEC, FUNBSEC = FUNBSEC,
                        FUNWEC = FUNWEC, FUNWSEC = FUNWSEC, FUNB = FUNB,
                        FUNW = FUNW)
      
      getcov <- getCov_JMH(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, 
                       H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS,
                       FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC,
                       FUNBSEC, FUNWEC, FUNWSEC,FUNB, FUNW)
      
      vcov <- getcov$vcov
      sebeta <- getcov$sebeta
      setau <- getcov$setau
      segamma1 <- getcov$segamma1
      segamma2 <- getcov$segamma2
      sealpha1 <- getcov$sealpha1
      sealpha2 <- getcov$sealpha2
      sevee1 <- getcov$sevee1
      sevee2 <- getcov$sevee2
      seSig <- getcov$seSig
      
      ### get loglike
      
      getloglike <- getLoglike.JMH(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, 
                               H01, H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, 
                               mdataS, xsmatrix, wsmatrix, method, initial.optimizer)
      
      names(beta) <- namesbeta
      names(tau) <- namestau
      names(gamma1) <- paste0(namesgamma1, "_1")
      names(gamma2) <- paste0(namesgamma1, "_2")
      
      FunCall_long <- longfmla
      FunCall_survival <- survfmla
      FunCall_longVar <- as.formula(paste("log(sigma^2)", varformula[2], sep = "~"))
      
      PropComp <- as.data.frame(table(cdata[, survival[2]]))
      
      ## return the joint modeling result
      mycall <- match.call()
      
      result <- list(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, 
                     H02, Sig, iter, convergence, vcov, sebeta, setau, segamma1,
                     segamma2, sealpha1, sealpha2, sevee1, sevee2, seSig, getloglike,
                     EFuntheta, CompetingRisk, quadpoint, rawydata, rawcdata, PropComp, 
                     FunCall_long, FunCall_longVar, FunCall_survival, random, method, mycall, epsilon)
      
      names(result) <- c("beta", "tau", "gamma1", "gamma2", "alpha1", "alpha2", "vee1",
                         "vee2", "H01", "H02", "Sig", "iter", "convergence", "vcov",
                         "sebeta", "setau", "segamma1", "segamma2", "sealpha1", "sealpha2", 
                         "sevee1", "sevee2", "seSig", "loglike", "EFuntheta",
                         "CompetingRisk", "quadpoint",
                         "ydata", "cdata", "PropEventType", "LongitudinalSubmodelmean",
                         "LongitudinalSubmodelvariance", "SurvivalSubmodel", "random", "method",
                         "call", "epsilon")
      
      class(result) <- "JMMLSM"
      
      return(result)
    }
  } else {
    repeat
    {
      iter <- iter + 1
      prebeta <- beta
      pretau <- tau
      pregamma1 <- gamma1
      prealpha1 <- alpha1
      prevee1 <- vee1
      preSig <- Sig
      preH01 <- H01
      if (print.para) {
        writeLines("iter is:")
        print(iter)
        writeLines("beta is:")
        print(beta)
        writeLines("tau is:")
        print(tau)
        writeLines("gamma1 is:")
        print(gamma1)
        writeLines("alpha1 is:")
        print(alpha1)
        writeLines("vee1 is:")
        print(vee1)
        writeLines("Sig is:")
        print(Sig)
      }
      
      if (method == "standard") {
        GetEfun <- GetESF.JMH(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, 
                          X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix)
      } else if (method == "adaptive") {
        GetEfun <- GetESFad.JMH(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, 
                            survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, initial.optimizer)
      } else {
        stop("Please choose one of the following methods for numerical integration in the E-step: standard, adaptive.")
      }
      
      
      GetMpara <- GetMSF.JMH(GetEfun, beta, tau, gamma1, alpha1, vee1, H01,
                         Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS)
      
      beta <- GetMpara$beta
      tau <- GetMpara$tau
      gamma1 <- GetMpara$gamma1
      alpha1 <- GetMpara$alpha1
      vee1 <- GetMpara$vee1
      Sig <- GetMpara$Sig
      H01 <- GetMpara$H01
      
      if((DiffSFrelative.JMH(beta, prebeta, tau, pretau, gamma1, pregamma1,
                 alpha1, prealpha1, vee1, prevee1, Sig, preSig, H01, preH01, epsilon) == 0) 
         || (iter == maxiter) || (!is.list(GetEfun)) || (!is.list(GetMpara))) {
        break
      }
      
    }
    
    if (iter == maxiter) {
      writeLines("program stops because of nonconvergence")
      convergence = 0
      result <- list(beta, tau, gamma1, alpha1, vee1, H01, Sig, iter, convergence)
      
      names(result) <- c("beta", "tau", "gamma1", "alpha1", "vee1", "H01", "Sig", 
                         "iter", "convergence")
      
      class(result) <- "JMMLSM"
      
      return(result)
    } else if (!is.list(GetEfun)) {
      writeLines("Something wrong in the E steps")
      beta <- NULL
      tau <- NULL
      gamma1 <- NULL
      alpha1 <- NULL
      vee1 <- NULL
      H01 <- NULL
      Sig <- NULL
      iter <- NULL
      result <- list(beta, tau, gamma1, alpha1, vee1, H01, Sig, iter)
      names(result) <- c("beta", "tau", "gamma1", "alpha1", "vee1", "H01", "Sig", "iter")
      class(result) <- "JMMLSM"
      return(result)
    } else if (!is.list(GetMpara)) {
      writeLines("Something wrong in the M steps")
      beta <- NULL
      tau <- NULL
      gamma1 <- NULL
      alpha1 <- NULL
      vee1 <- NULL
      H01 <- NULL
      Sig <- NULL
      iter <- NULL
      result <- list(beta, tau, gamma1, alpha1, vee1, H01, Sig, iter)
      names(result) <- c("beta", "tau", "gamma1", "alpha1", "vee1", "H01", "Sig", "iter")
      class(result) <- "JMMLSM"
      return(result)
    } else {
      
      convergence = 1
      
      if (method == "standard") {
        GetEfun <- GetESF.JMH(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, 
                          X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix)
      } else {
        GetEfun <- GetESFad.JMH(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, 
                            survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, initial.optimizer)
      }
      
      FUNENW <- as.vector(GetEfun$FUNENW)
      FUNBENW <- as.matrix(GetEfun$FUNBENW)
      FUNBS <- as.matrix(GetEfun$FUNBS)
      FUNBW <- as.matrix(GetEfun$FUNBW)
      FUNWS <- as.vector(GetEfun$FUNWS)
      FUNBSENW <- as.matrix(GetEfun$FUNBSENW) 
      FUNEC <- as.matrix(GetEfun$FUNEC)
      FUNBEC <- as.matrix(GetEfun$FUNBEC)
      FUNBSEC <- as.matrix(GetEfun$FUNBSEC)
      FUNWEC <- as.matrix(GetEfun$FUNWEC)
      FUNWSEC <- as.matrix(GetEfun$FUNWSEC)
      FUNB <- as.matrix(GetEfun$FUNB)
      FUNW <- as.vector(GetEfun$FUNW)
      
      getcov <- getCovSF_JMH(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y,
                         X2, survtime, cmprsk, mdata, mdataS,
                         FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC,
                         FUNBSEC, FUNWEC, FUNWSEC,FUNB, FUNW)
      
      vcov <- getcov$vcov
      seSig <- getcov$seSig
      sebeta <- vector()
      setau <- vector()
      segamma1 <- vector()
      sealpha1 <- vector()
      for (i in 1:length(beta)) sebeta[i] <- sqrt(vcov[i, i])
      for (i in 1:length(tau)) setau[i] <- sqrt(vcov[length(beta)+i, length(beta)+i])
      for (i in 1:length(gamma1)) segamma1[i] <- sqrt(vcov[length(beta)+length(tau)+i, 
                                                           length(beta)+length(tau)+i])
      for (i in 1:length(alpha1)) sealpha1[i] <- sqrt(vcov[length(beta)+length(tau)+length(gamma1)+i, 
                                                           length(beta)+length(tau)+length(gamma1)+i])
      sevee1 <- sqrt(vcov[length(beta)+length(tau)+length(gamma1)+length(alpha1)+1, 
                          length(beta)+length(tau)+length(gamma1)+length(alpha1)+1])
      
      
      getloglike <- getLoglikeSF.JMH(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, 
                                 X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, method, initial.optimizer)
      
      
      names(beta) <- namesbeta
      names(tau) <- namestau
      names(gamma1) <- paste0(namesgamma1, "_1")
      
      FunCall_long <- longfmla
      FunCall_survival <- survfmla
      FunCall_longVar <- as.formula(paste("log(sigma^2)", varformula[2], sep = "~"))
      
      PropComp <- as.data.frame(table(cdata[, survival[2]]))
      
      ## return the joint modeling result
      mycall <- match.call()
      
      result <- list(beta, tau, gamma1, alpha1, vee1, H01, Sig, iter, convergence, 
                     vcov, sebeta, setau, segamma1, sealpha1, sevee1, seSig, getloglike,
                     CompetingRisk, quadpoint, rawydata, rawcdata, PropComp, 
                     FunCall_long, FunCall_longVar, FunCall_survival, random, mycall, method, epsilon)
      
      names(result) <- c("beta", "tau", "gamma1", "alpha1", "vee1", "H01", "Sig", 
                         "iter", "convergence", "vcov", "sebeta", "setau", "segamma1", 
                         "sealpha1", "sevee1", "seSig", "loglike", "CompetingRisk", "quadpoint",
                         "ydata", "cdata", "PropEventType", "LongitudinalSubmodelmean",
                         "LongitudinalSubmodelvariance", "SurvivalSubmodel", "random",
                         "call", "method", "epsilon")
      
      class(result) <- "JMMLSM"
      
      return(result)
    }
    
    
    
    
  }
  
  
  
  
}