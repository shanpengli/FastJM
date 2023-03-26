##' Joint modeling of longitudinal continuous data and competing risks
##' @title Joint modeling of longitudinal continuous data and competing risks
##' @name jmcs
##' @param ydata a longitudinal data frame in long format.
##' @param cdata a survival data frame with competing risks or single failure.
##' Each subject has one data entry.
##' @param long.formula a formula object with the response variable and fixed effects covariates
##' to be included in the longitudinal sub-model.
##' @param random a one-sided formula object describing the random effects part of the longitudinal sub-model.
##' For example, fitting a random intercept model takes the form \code{ ~ 1|ID}.
##' Alternatively. Fitting a random intercept and slope model takes the form \code{~ x1 + ... + xn|ID}.
##' @param surv.formula a formula object with the survival time, event indicator, and the covariates
##' to be included in the survival sub-model.
##' @param REML a logic object that indicates the use of REML estimator. Default is TRUE.
##' @param quadpoint the number of pseudo-adaptive Gauss-Hermite quadrature points.
##' to be chosen for numerical integration. Default is 6 which produces stable estimates in most dataframes.
##' @param maxiter the maximum number of iterations of the EM algorithm that the function will perform. Default is 10000.
##' @param print.para Print detailed information of each iteration. Default is FALSE, i.e., not to print the iteration details.
##' @param survinitial Fit a Cox model to obtain initial values of the parameter estimates. Default is TRUE.
##' @param tol Tolerance parameter. Default is 0.0001.
##' @param method Method for proceeding numerical integration in the E-step. Default is pseudo-adaptive. 
##' @param opt Optimization method to fit a linear mixed effects model, either \code{nlminb} (default) or \code{optim}.
##' @return  Object of class \code{jmcs} with elements
##' \item{beta}{the vector of fixed effects for the linear mixed effects model.} 
##' \item{gamma1}{the vector of fixed effects for type 1 failure for the survival model.}
##' \item{gamma2}{the vector of fixed effects for type 2 failure for the survival model. 
##' Valid only if \code{CompetingRisk = TRUE}.}
##' \item{nu1}{the vector of association parameter(s) for type 1 failure.}
##' \item{nu2}{the vector of association parameter(s) for type 2 failure. Valid only if \code{CompetingRisk = TRUE}.}
##' \item{H01}{the matrix that collects baseline hazards evaluated at each uncensored event time for type 1 failure. 
##' The first column denotes uncensored event times, the second column the number of events, and the third columns 
##' the hazards obtained by Breslow estimator.}
##' \item{H02}{the matrix that collects baseline hazards evaluated at each uncensored event time for type 2 failure. 
##' The data structure is the same as \code{H01}. Valid only if \code{CompetingRisk = TRUE}.}
##' \item{Sig}{the variance-covariance matrix of the random effects.}
##' \item{sigma}{the variance of the measurement error for the linear mixed effects model.}
##' \item{iter}{the total number of iterations until convergence.}
##' \item{convergence}{convergence identifier: 1 corresponds to successful convergence, 
##' whereas 0 to a problem (i.e., when 0, usually more iterations are required).}
##' \item{vcov}{the variance-covariance matrix of all the fixed effects for both models.}
##' \item{sebeta}{the standard error of \code{beta}.}
##' \item{segamma1}{the standard error of \code{gamma1}.}
##' \item{segamma2}{the standard error of \code{gamma2}. 
##' Valid only if \code{CompetingRisk = TRUE}.}
##' \item{senu1}{the standard error of \code{nu1}.}
##' \item{senu2}{the standard error of \code{nu2}. Valid only if \code{CompetingRisk = TRUE}.}
##' \item{seSig}{the vector of standard errors of covariance of random effects.}
##' \item{sesigma}{the standard error of variance of measurement error for the linear mixed effects model.}
##' \item{loglike}{the log-likelihood value.}
##' \item{fitted}{a list with the fitted values:
##'   \describe{
##'   \item{resid}{the vector of estimated residuals for the linear mixed effects model.} 
##'   \item{fitted}{the vector of fitted values for the linear mixed effects model.}
##'   \item{fittedmar}{the vector of marginal fitted values for the linear mixed effects model.}
##'   \item{residmar}{the vector of estimated marginal residuals for the linear mixed effects model.}
##'   }
##' }
##' \item{fittedSurv}{the estimated survival rate evaluated at each uncensored event time.}
##' \item{FUNB}{the estimated random effects for each subject.}
##' \item{CompetingRisk}{logical value; TRUE if a competing event are accounted for.}
##' \item{quadpoint}{the number of Gauss Hermite quadrature points used for numerical integration.}
##' \item{ydata}{the input longitudinal dataset for fitting a joint model.
##' It has been re-ordered in accordance with descending observation times in \code{cdata}.}
##' \item{cdata}{the input survival dataset for fitting a joint model.
##' It has been re-ordered in accordance with descending observation times.}
##' \item{PropEventType}{a frequency table of number of events.}
##' \item{LongitudinalSubmodel}{the component of the \code{long.formula}.}
##' \item{SurvivalSubmodel}{the component of the \code{surv.formula}.}
##' \item{random}{the component of the \code{random}.}
##' \item{call}{the matched call.}
##' \item{Quad.method}{the quadrature rule used for integration. 
##' If pseudo-adaptive quadrature rule is used, then return \code{pseudo-adaptive}. 
##' Otherwise return \code{standard}.}
##' \item{id}{the grouping vector for the longitudinal outcome.}
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{ranef}, \link{fixef}, \link{fitted.jmcs}, 
##' \link{residuals.jmcs}, \link{survfitjmcs}, \link{plot.jmcs},
##' \link{vcov.jmcs}}
##' @examples 
##' 
##' require(FastJM)
##' # Load a simulated longitudinal dataset
##' data(ydata)
##' # Load a simulated survival dataset with two competing events
##' data(cdata)
##' \donttest{
##' # Fit a joint model
##' fit <- jmcs(ydata = ydata, cdata = cdata, 
##'             long.formula = response ~ time + gender + x1 + race, 
##'             surv.formula = Surv(surv, failure_type) ~ x1 + gender + x2 + race, 
##'             random =  ~ time| ID)
##' fit
##' # Extract the parameter estimates of longitudinal sub-model fixed effects
##' fixef(fit, process = "Longitudinal")
##' # Extract the parameter estimates of survival sub-model fixed effects
##' fixef(fit, process = "Event")
##' # Obtain the random effects estimates for first 6 subjects 
##' head(ranef(fit))
##' # Obtain the variance-covariance matrix of all parameter estimates 
##' vcov(fit)
##' # Obtain the result summaries of the joint model fit
##' summary(fit, process = "Longitudinal")
##' summary(fit, process = "Event")
##' # Prediction of cumulative incidence for competing risks data
##' # Predict the conditional probabilities for two patients who are alive (censored)
##' ND <- ydata[ydata$ID %in% c(419, 218), ]
##' ID <- unique(ND$ID)
##' NDc <- cdata[cdata$ID  %in% ID, ]
##' survfit <- survfitjmcs(fit, 
##'                        ynewdata = ND, 
##'                        cnewdata = NDc, 
##'                        u = seq(3, 4.8, by = 0.2), 
##'                        method = "GH",
##'                        obs.time = "time")
##' survfit
##' PE <- PEjmcs(fit, seed = 100, landmark.time = 3, horizon.time = c(3.6, 4, 4.4), 
##'              obs.time = "time", method = "GH", 
##'              quadpoint = NULL, maxiter = 1000, n.cv = 3, 
##'              survinitial = TRUE)
##' Brier <- summary(PE, error = "Brier")
##' Brier
##' 
##' MAEQ <- MAEQjmcs(fit, seed = 100, landmark.time = 3, horizon.time = c(3.6, 4, 4.4), 
##'                  obs.time = "time", method = "GH", 
##'                  quadpoint = NULL, maxiter = 1000, n.cv = 3, 
##'                  survinitial = TRUE)
##' APE <- summary(MAEQ, digits = 3)
##' APE
##' }
##' 
##' @export
##'
##'

jmcs <- function(ydata, cdata, long.formula, random = NULL, surv.formula, REML = TRUE,
                 quadpoint = NULL, maxiter = 10000, print.para = FALSE, survinitial = TRUE, tol = 0.0001, 
                 method = "pseudo-adaptive", opt = "nlminb")
{
  
  if (method == "pseudo-adaptive" & is.null(quadpoint)) {
    quadpoint <- 6
  }
  
  if (method == "standard" & is.null(quadpoint)) {
    quadpoint <- 20
  }
  
  if (!(method %in% c("standard", "pseudo-adaptive")))
    stop("Please choose one of the following metho for numerical integration: standard, pseudo-adaptive. See help() for details of the methods.")
  
  if (!inherits(long.formula, "formula") || length(long.formula) != 3) {
    stop("\nLinear mixed effects model must be a formula of the form \"resp ~ pred\"")
  }
  
  if (!inherits(surv.formula, "formula") || length(surv.formula) != 3) {
    stop("\nCox proportional hazards model must be a formula of the form \"Surv(.,.) ~ pred\"")
  }
  
  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
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

  getinit <- Getinit(cdata = cdata, ydata = ydata, long.formula = long.formula,
                     surv.formula = surv.formula,
                     model = model, ID = ID, RE = RE, survinitial = survinitial, 
                     REML = REML, random = random, opt = opt)
  
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
    gamma1 <- getinit$gamma1
    namesgamma1 <- names(gamma1)
    gamma2 <- getinit$gamma2
    alpha1 <- getinit$alpha1
    alpha2 <- getinit$alpha2
    Sig <- getinit$Sig
    p1a <- ncol(Sig)
    
    CompetingRisk <- TRUE
  } else {
    ## initialize parameters
    getriskset <- Getriskset(cdata = cdata, surv.formula = surv.formula)
    
    ## number of distinct survival time
    H01 <- as.matrix(getriskset$tablerisk1)
    
    ## initialize parameters
    beta <- getinit$beta
    namesbeta <- names(beta)
    gamma1 <- getinit$gamma1
    namesgamma1 <- names(gamma1)
    alpha1 <- getinit$alpha1
    Sig <- getinit$Sig
    p1a <- ncol(Sig)
    
    CompetingRisk <- FALSE
  }
  
  if (p1a > 3 | p1a < 1) {
    stop("The current package cannot handle this dimension of random effects. Please re-consider your model.")
  }
  
  getGH <- GetGHmatrix(quadpoint = quadpoint, Sig = Sig)
  
  xsmatrix <- getGH$xsmatrix
  wsmatrix <- getGH$wsmatrix
  
  iter=0
  
  Z <- getinit$Z
  X1 <- getinit$X1
  Y <- getinit$Y
  X2 <- getinit$X2
  sigma <- getinit$sigma
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
  
  if (method == "pseudo-adaptive") {

    getBayesEst <- getBayes(beta, Sig, sigma, Z, X1, Y, mdata, mdataS)
    Posbi <- getBayesEst$Posbi
    Poscov <- getBayesEst$Poscov
  } else {
    Posbi <- 0
    Poscov <- 0
  }

  
  if (CompetingRisk == TRUE) {
    repeat
    {
      iter <- iter + 1
      prebeta <- beta
      pregamma1 <- gamma1
      pregamma2 <- gamma2
      prealpha1 <- alpha1
      prealpha2 <- alpha2
      preSig <- Sig
      presigma <- sigma
      preH01 <- H01
      preH02 <- H02
      if (print.para) {
        writeLines("iter is:")
        print(iter)
        writeLines("beta is:")
        print(beta)
        writeLines("gamma1 is:")
        print(gamma1)
        writeLines("gamma2 is:")
        print(gamma2)
        writeLines("nu1 is:")
        print(alpha1)
        writeLines("nu2 is:")
        print(alpha2)
        writeLines("Sig is:")
        print(Sig)
        writeLines("Error variance is:")
        print(sigma)
      }
      
      GetEfun <- GetE(beta, gamma1, gamma2, alpha1, alpha2, H01, H02, 
                      Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, method, Posbi, Poscov)

      GetMpara <- GetM(GetEfun, beta, gamma1, gamma2, alpha1, alpha2, H01, H02, 
                       Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS)

      
      beta <- GetMpara$beta
      sigma <- GetMpara$sigma
      gamma1 <- GetMpara$gamma1
      gamma2 <- GetMpara$gamma2
      alpha1 <- GetMpara$alpha1
      alpha2 <- GetMpara$alpha2
      Sig <- GetMpara$Sig
      sigma <- GetMpara$sigma
      H01 <- GetMpara$H01
      H02 <- GetMpara$H02
      
      if((Diff(beta, prebeta, sigma, presigma, gamma1, pregamma1, gamma2, pregamma2,
               alpha1, prealpha1, alpha2, prealpha2,
               Sig, preSig, H01, preH01, H02, preH02, tol) == 0) || (iter == maxiter) || (!is.list(GetEfun))
         || (!is.list(GetMpara))) {
        break
      }
      
    }
    
    if (iter == maxiter) {
      writeLines("program stops because of nonconvergence")
      convergence = 0
      result <- list(beta, gamma1, gamma2, alpha1, alpha2, H01, 
                     H02, Sig, sigma, iter, convergence)
      
      names(result) <- c("beta", "gamma1", "gamma2", "nu1", "nu2", 
                         "H01", "H02", "Sig", "sigma", "iter", "convergence")
      
      return(result)
    } else if (!is.list(GetEfun)) {
      writeLines("Something wrong in the E steps")
      beta <- NULL
      gamma1 <- NULL
      gamma2 <- NULL
      alpha1 <- NULL
      alpha2 <- NULL
      H01 <- NULL
      H02 <- NULL
      Sig <- NULL
      sigma <- NULL
      iter <- NULL
      result <- list(beta, gamma1, gamma2, alpha1, alpha2, H01, 
                     H02, Sig, sigma, iter)
      
      names(result) <- c("beta", "gamma1", "gamma2", "nu1", "nu2", 
                         "H01", "H02", "Sig", "sigma", "iter")
      
      return(result)
    } else if (!is.list(GetMpara)) {
      writeLines("Something wrong in the M steps")
      beta <- NULL
      gamma1 <- NULL
      gamma2 <- NULL
      alpha1 <- NULL
      alpha2 <- NULL
      H01 <- NULL
      H02 <- NULL
      Sig <- NULL
      sigma <- NULL
      iter <- NULL
      result <- list(beta, gamma1, gamma2, alpha1, alpha2, H01, 
                     H02, Sig, sigma, iter)
      names(result) <- c("beta", "gamma1", "gamma2", "nu1", "nu2", 
                         "H01", "H02", "Sig", "sigma",  "iter")
      
      return(result)
    } else {
      
      convergence = 1
      
      GetEfun <- GetE(beta, gamma1, gamma2, alpha1, alpha2, H01, H02, 
                      Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, method, Posbi, Poscov)
      
      FUNBS <- as.matrix(GetEfun$FUNBS)
      FUNEC <- as.matrix(GetEfun$FUNEC)
      FUNBEC <- as.matrix(GetEfun$FUNBEC)
      FUNBSEC <- as.matrix(GetEfun$FUNBSEC)
      FUNB <- as.matrix(GetEfun$FUNB)
      
      getcov <- getCov(beta, gamma1, gamma2, alpha1, alpha2, H01, 
                       H02, Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS,
                       FUNBS, FUNEC, FUNBEC, FUNBSEC, FUNB)
      
      vcov <- getcov$vcov
      sebeta <- getcov$sebeta
      segamma1 <- getcov$segamma1
      segamma2 <- getcov$segamma2
      sealpha1 <- getcov$sealpha1
      sealpha2 <- getcov$sealpha2
      seSig <- getcov$seSig
      sesigma <- getcov$sesigma
      
      getloglike <- getLoglike(beta, gamma1, gamma2, alpha1, alpha2, H01, H02, 
                               Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS, 
                               xsmatrix, wsmatrix, method, Posbi, Poscov)
      
      getfitted <- getfitted(beta, Z, X1, Y, mdata, mdataS, FUNB)
      
      CH01 <- data.frame(H01[, 1], cumsum(H01[, 3]), NA)
      colnames(CH01) <- c("Time", "CH01", "CH02")
      CH02 <- data.frame(H02[, 1], NA, cumsum(H02[, 3]))
      colnames(CH02) <- c("Time", "CH01", "CH02")
      CH012 <- rbind(CH01, CH02)
      CH012 <- CH012[order(CH012$Time), ]
      if (is.na(CH012$CH01[1])) CH012$CH01[1] <- 0
      if (is.na(CH012$CH02[1])) CH012$CH02[1] <- 0
      CH012[is.na(CH012)] <- 0 
      CH012 <- as.matrix(CH012)
      getfittedSurv <- getfittedSurv(gamma1, gamma2, X2, CH012, alpha1, alpha2, FUNB)
      
      survival <- all.vars(surv.formula)
      id <- ydata[, ID]
      
      names(beta) <- namesbeta
      
      names(gamma1) <- paste0(namesgamma1, "_1")
      names(gamma2) <- paste0(namesgamma1, "_2")
      
      PropComp <- as.data.frame(table(cdata[, survival[2]]))
      
      FunCall_long <- longfmla
      
      FunCall_survival <- survfmla
      
      ## return the joint modelling result
      mycall <- match.call()
      
      result <- list(beta, gamma1, gamma2, alpha1, alpha2, H01, 
                     H02, Sig, sigma, iter, convergence, vcov, sebeta, segamma1,
                     segamma2, sealpha1, sealpha2, seSig, sesigma, getloglike, 
                     getfitted, getfittedSurv, FUNB, CompetingRisk,
                     quadpoint, ydata, cdata, PropComp, FunCall_long,
                     FunCall_survival, random, mycall, method, id)
      
      names(result) <- c("beta", "gamma1", "gamma2", "nu1", "nu2", "H01", "H02", "Sig", 
                         "sigma", "iter", "convergence", "vcov",
                         "sebeta", "segamma1", "segamma2", "senu1", "senu2", 
                          "seSig", "sesigma", "loglike", "fitted", "fittedSurv", 
                         "FUNB", "CompetingRisk", "quadpoint",
                         "ydata", "cdata", "PropEventType", "LongitudinalSubmodel",
                          "SurvivalSubmodel", "random", "call", "Quad.method", "id")
      
      class(result) <- "jmcs"
      
      return(result)
    }
  } else {
    repeat
    {
      iter <- iter + 1
      prebeta <- beta
      pregamma1 <- gamma1
      prealpha1 <- alpha1
      preSig <- Sig
      presigma <- sigma
      preH01 <- H01
      if (print.para) {
        writeLines("iter is:")
        print(iter)
        writeLines("beta is:")
        print(beta)
        writeLines("gamma1 is:")
        print(gamma1)
        writeLines("nu1 is:")
        print(alpha1)
        writeLines("Sig is:")
        print(Sig)
        writeLines("sigma is:")
        print(sigma)
      }
      
      GetEfun <- GetESF(beta, gamma1, alpha1, H01, Sig, sigma, Z, X1, Y, 
                        X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, method, Posbi, Poscov)
      
      
      GetMpara <- GetMSF(GetEfun, beta, gamma1, alpha1, H01,
                         Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS)
      
      beta <- GetMpara$beta
      gamma1 <- GetMpara$gamma1
      alpha1 <- GetMpara$alpha1
      Sig <- GetMpara$Sig
      sigma <- GetMpara$sigma
      H01 <- GetMpara$H01
      
      if((DiffSF(beta, prebeta, sigma, presigma, gamma1, pregamma1,
                 alpha1, prealpha1, Sig, preSig, H01, preH01, tol) == 0) 
         || (iter == maxiter) || (!is.list(GetEfun)) || (!is.list(GetMpara))) {
        break
      }
      
    }
    
    if (iter == maxiter) {
      writeLines("program stops because of nonconvergence")
      convergence = 0
      result <- list(beta, gamma1, alpha1, H01, Sig, sigma, iter, convergence)
      
      names(result) <- c("beta", "gamma1", "nu1", "H01", "Sig", "sigma", 
                         "iter", "convergence")
      
      return(result)
    } else if (!is.list(GetEfun)) {
      writeLines("Something wrong in the E steps")
      beta <- NULL
      gamma1 <- NULL
      alpha1 <- NULL
      H01 <- NULL
      Sig <- NULL
      sigma <- NULL
      iter <- NULL
      result <- list(beta, gamma1, alpha1, H01, Sig, sigma, iter)
      names(result) <- c("beta", "gamma1", "nu1", "H01", "Sig", "sigma", "iter")
      return(result)
    } else if (!is.list(GetMpara)) {
      writeLines("Something wrong in the M steps")
      beta <- NULL
      gamma1 <- NULL
      alpha1 <- NULL
      H01 <- NULL
      Sig <- NULL
      sigma <- NULL
      iter <- NULL
      result <- list(beta, gamma1, alpha1, H01, Sig, sigma, iter)
      names(result) <- c("beta", "gamma1", "nu1", "H01", "Sig", "sigma", "iter")
      return(result)
    } else {
      
      convergence = 1
      
      GetEfun <- GetESF(beta, gamma1, alpha1, H01, Sig, sigma, Z, X1, Y, 
                        X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, method, Posbi, Poscov)
      
      FUNBS <- as.matrix(GetEfun$FUNBS)
      FUNEC <- as.matrix(GetEfun$FUNEC)
      FUNBEC <- as.matrix(GetEfun$FUNBEC)
      FUNBSEC <- as.matrix(GetEfun$FUNBSEC)
      FUNB <- as.matrix(GetEfun$FUNB)
      
      getcov <- getCovSF(beta, gamma1, alpha1, H01, Sig, sigma, Z, X1, Y, 
                         X2, survtime, cmprsk, mdata, mdataS,
                         FUNBS, FUNEC, FUNBEC, FUNBSEC, FUNB)
      
      vcov <- getcov$vcov
      sebeta <- getcov$sebeta
      segamma1 <- getcov$segamma1
      sealpha1 <- getcov$sealpha1
      sesigma <- getcov$sesigma
      seSig <- getcov$seSig
      
      getloglike <- getLoglikeSF(beta, gamma1, alpha1, H01,
                               Sig, sigma, Z, X1, Y, X2, 
                               survtime, cmprsk, mdata, mdataS, 
                               xsmatrix, wsmatrix, method, Posbi, Poscov)
      
      getfitted <- getfitted(beta, Z, X1, Y, mdata, mdataS, FUNB)
      
      CH01 <- data.frame(H01, cumsum(H01[, 3]))
      CH01 <- as.matrix(CH01)
      getfittedSurv <- getfittedSurvSF(gamma1, X2, CH01, alpha1, FUNB)
      
      survival <- all.vars(surv.formula)
      id <- ydata[, ID]
      
      names(beta) <- namesbeta 
      names(gamma1) <- paste0(namesgamma1, "_1")
      # names(gamma1) <- paste0(survival[-(1:2)], "_1")
      
      PropComp <- as.data.frame(table(cdata[, survival[2]]))
      
      survival <- all.vars(survfmla)
      
      FunCall_long <- longfmla
      
      # SurvOut <- paste0("Surv(", survival[1], ",", survival[2], ")")
      # SurvX <- paste0(survival[-(1:2)], collapse = "+")
      # FunCall_survival <- as.formula(paste(SurvOut, SurvX, sep = "~"))
      # 
      FunCall_survival <- survfmla
      
      ## return the joint modelling result
      mycall <- match.call()
      
      result <- list(beta, gamma1, alpha1, H01, Sig, sigma, iter, convergence, 
                     vcov, sebeta, segamma1, sealpha1, seSig, sesigma, getloglike, 
                     getfitted, getfittedSurv, FUNB, CompetingRisk,
                     quadpoint, ydata, cdata, PropComp, FunCall_long, FunCall_survival, 
                     random, mycall, method, id)
      
      names(result) <- c("beta", "gamma1", "nu1", "H01", "Sig", "sigma",
                         "iter", "convergence", "vcov", "sebeta", "segamma1", 
                         "senu1", "seSig", "sesigma", "loglike", "fitted", "fittedSurv", 
                         "FUNB", "CompetingRisk", "quadpoint",
                         "ydata", "cdata", "PropEventType", "LongitudinalSubmodel",
                         "SurvivalSubmodel", "random", "call", "Quad.method", "id")
      
      class(result) <- "jmcs"
      
      return(result)
    }
    
    
    
    
  }
  
  
  
  
  
  
  
}