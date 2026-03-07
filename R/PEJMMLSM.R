##' @title A metric of prediction accuracy of joint model by comparing the predicted risk
##' with the counting process.
##' @name PEJMMLSM
##' @aliases PEJMMLSM
##' @param seed a numeric value of seed to be specified for cross validation.
##' @param object object of class 'JMMLSM'.
##' @param landmark.time a numeric value of time for which dynamic prediction starts..
##' @param horizon.time a numeric vector of future times for which predicted probabilities are to be computed.
##' @param obs.time a character string of specifying a longitudinal time variable.
##' @param method estimation method for predicted probabilities. If \code{Laplace}, then the empirical empirical
##' estimates of random effects is used. If \code{GH}, then the pseudo-adaptive Gauss-Hermite quadrature is used.
##' @param quadpoint the number of pseudo-adaptive Gauss-Hermite quadrature points if \code{method = "GH"}.
##' @param maxiter the maximum number of iterations of the EM algorithm that the 
##' function will perform. Default is 10000.
##' @param n.cv number of folds for cross validation. Default is 3.
##' @param survinitial Fit a Cox model to obtain initial values of the parameter estimates. Default is TRUE.
##' @param opt Optimization method to fit a linear mixed effects model, either nlminb (default) or optim.
##' @param initial.para Initial guess of parameters for cross validation. Default is FALSE.
##' @param LOCF a logical value to indicate whether the last-observation-carried-forward approach applies to prediction. 
##' If \code{TRUE}, then \code{LOCFcovariate} and \code{clongdata} must be specified to indicate 
##' which time-dependent survival covariates are included for dynamic prediction. Default is FALSE.
##' @param LOCFcovariate a vector of string with time-dependent survival covariates if \code{LOCF = TRUE}. Default is NULL.
##' @param clongdata a long format data frame where time-dependent survival covariates are incorporated. Default is NULL.
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{JMMLSM}, \link{survfitJMMLSM}}
##' @export
##' 

PEJMMLSM <- function(seed = 100, object, landmark.time = NULL, horizon.time = NULL, 
                  obs.time = NULL, method = c("Laplace", "GH"), 
                  quadpoint = NULL, maxiter = 1000, n.cv = 3, survinitial = TRUE, 
                  opt = "nlminb", initial.para = FALSE, 
                  LOCF = FALSE, LOCFcovariate = NULL, clongdata = NULL, ...) {
  
  if (!inherits(object, "JMMLSM"))
    stop("Use only with 'JMMLSM' xs.\n")
  if (is.null(landmark.time)) 
    stop("Please specify the landmark.time for dynamic prediction.")   
  if (!method %in% c("Laplace", "GH"))
    stop("Please specify a method for probability approximation: Laplace or GH.")
  if (!is.vector(horizon.time)) 
    stop("horizon.time must be vector typed.")
  if (is.null(quadpoint)) {
    quadpoint <- object$quadpoint
  }
  if (is.null(obs.time)) {
    stop("Please specify a vector that represents the time variable from ydatanew.")
  } else {
    if (!obs.time %in% colnames(object$ydata)) {
      stop(paste0(obs.time, " is not found in ynewdata."))
    }
  }
  CompetingRisk <- object$CompetingRisk
  set.seed(seed)
  cdata <- object$cdata
  ydata <- object$ydata
  long.formula <- object$LongitudinalSubmodelmean
  surv.formula <- object$SurvivalSubmodel
  surv.var <- all.vars(surv.formula)
  New.surv.formula.out <- paste0("survival::Surv(", surv.var[1], ",", 
                                 surv.var[2], "==0)")
  New.surv.formula <- as.formula(paste(New.surv.formula.out, 1, sep = "~"))
  variance.formula <- as.formula(paste("", object$LongitudinalSubmodelvariance[3], sep = "~"))
  random <- all.vars(object$random) 
  ID <- random[length(random)]
  
  if (initial.para) {
    initial.para <- list(beta = object$beta,
                         tau = object$tau, 
                         gamma1 = object$gamma1,
                         gamma2 = object$gamma2,
                         alpha1 = object$alpha1,
                         alpha2 = object$alpha2,
                         vee1 = object$vee1,
                         vee2 = object$vee2,
                         Sig = object$Sig)
  } else {
    initial.para <- NULL
  }
  
  folds <- caret::groupKFold(c(1:nrow(cdata)), k = n.cv)
  Brier.cv <- list()
  MAE.cv <- list()
  for (t in 1:n.cv) {
    
    train.cdata <- cdata[folds[[t]], ]
    train.ydata <- ydata[ydata[, ID] %in% train.cdata[, ID], ]
    
    fit <- try(JMMLSM(cdata = train.cdata, ydata = train.ydata, 
                      long.formula = long.formula,
                      surv.formula = surv.formula,
                      variance.formula = variance.formula, 
                      quadpoint = quadpoint, random = object$random, 
                      survinitial = survinitial, maxiter = maxiter, opt = opt, 
                      epsilon = object$epsilon,
                      initial.para = initial.para), silent = TRUE)
    
    if ('try-error' %in% class(fit)) {
      writeLines(paste0("Error occured in the ", t, " th training!"))
      Brier.cv[[t]] <- NULL
      MAE.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      Brier.cv[[t]] <- NULL
      MAE.cv[[t]] <- NULL
    } else {
      
      val.cdata <- cdata[-folds[[t]], ]
      val.ydata <- ydata[ydata[, ID] %in% val.cdata[, ID], ]
      ## fit a Kalplan-Meier estimator
      fitKM <- survival::survfit(New.surv.formula, data = val.cdata)
      
      val.cdata <- val.cdata[val.cdata[, surv.var[1]] > landmark.time, ]
      val.ydata <- val.ydata[val.ydata[, ID] %in% val.cdata[, ID], ]
      val.ydata <- val.ydata[val.ydata[, obs.time] <= landmark.time, ]
      NewyID <- unique(val.ydata[, ID])
      val.cdata <- val.cdata[val.cdata[, ID] %in% NewyID, ]
      
      if (LOCF) {
        val.clongdata <- clongdata[clongdata[, ID] %in% val.cdata[, ID], ]
      } else {
        val.clongdata <- NULL
      }
      
      survfit <- try(survfitJMMLSM(fit, ynewdata = val.ydata, cnewdata = val.cdata, 
                                   u = horizon.time, method = method, 
                                   Last.time = landmark.time,
                                   obs.time = obs.time, quadpoint = quadpoint,
                                   LOCF = LOCF, LOCFcovariate = LOCFcovariate, clongdata = val.clongdata), silent = TRUE)
      
      if ('try-error' %in% class(survfit)) {
        writeLines(paste0("Error occured in the ", t, " th validation!"))
        Brier.cv[[t]] <- NULL
        MAE.cv[[t]] <- NULL
      } else {
        if (CompetingRisk) {
          CIF <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 3))
          colnames(CIF) <- c("ID", "CIF1", "CIF2")
          CIF$ID <- val.cdata[, ID]
          Gs <- summary(fitKM, times = landmark.time)$surv
          mean.Brier <- matrix(NA, nrow = length(horizon.time), ncol = 2)
          mean.MAE <- matrix(NA, nrow = length(horizon.time), ncol = 2)
          for (j in 1:length(horizon.time)) {
            fitKM.horizon <- try(summary(fitKM, times = horizon.time[j]), silent = TRUE)
            if ('try-error' %in% class(fitKM.horizon)) {
              mean.Brier[j, 1] <- NA
              mean.Brier[j, 2] <- NA
              mean.MAE[j, 1] <- NA
              mean.MAE[j, 2] <- NA
            } else {
              Gu <- fitKM.horizon$surv
              ## true counting process
              N1 <- vector()
              N2 <- vector()
              Gt <- vector()
              W.IPCW <- vector()
              for (i in 1:nrow(CIF)) {
                if (val.cdata[i, surv.var[1]] <= horizon.time[j] && val.cdata[i, surv.var[2]] == 1) {
                  N1[i] <- 1
                } else {
                  N1[i] <- 0
                }
                
                if (val.cdata[i, surv.var[1]] <= horizon.time[j] && val.cdata[i, surv.var[2]] == 2) {
                  N2[i] <- 1
                } else {
                  N2[i] <- 0
                }
                
                if (val.cdata[i, surv.var[1]] <= horizon.time[j] && val.cdata[i, surv.var[2]] != 0) {
                  Gt[i] <- summary(fitKM, times = val.cdata[i, surv.var[1]])$surv
                } else {
                  Gt[i] <- NA
                }
                
                if (val.cdata[i, surv.var[1]] > horizon.time[j]) {
                  W.IPCW[i] <- 1/(Gu/Gs)
                } else if (val.cdata[i, surv.var[1]] <= horizon.time[j] && val.cdata[i, surv.var[2]] != 0) {
                  W.IPCW[i] <- 1/(Gt[i]/Gs)
                } else {
                  W.IPCW[i] <- NA
                }
              }
              ## extract estimated CIF
              for (k in 1:nrow(CIF)) {
                CIF[k, 2] <- survfit$Pred[[k]][j, 2]
                CIF[k, 3] <- survfit$Pred[[k]][j, 3]
              }
              
              RAWData.Brier <- data.frame(CIF, N1, N2, W.IPCW)
              colnames(RAWData.Brier)[1:3] <- c("ID", "CIF1", "CIF2")
              RAWData.Brier$Brier1 <- RAWData.Brier$W.IPCW*
                abs(RAWData.Brier$CIF1 - RAWData.Brier$N1)^2
              RAWData.Brier$Brier2 <- RAWData.Brier$W.IPCW*
                abs(RAWData.Brier$CIF2 - RAWData.Brier$N2)^2
              RAWData.Brier$MAE1 <- RAWData.Brier$W.IPCW*
                abs(RAWData.Brier$CIF1 - RAWData.Brier$N1)^1
              RAWData.Brier$MAE2 <- RAWData.Brier$W.IPCW*
                abs(RAWData.Brier$CIF2 - RAWData.Brier$N2)^1
              mean.MAE1  <- sum(RAWData.Brier$MAE1, na.rm = TRUE)/nrow(RAWData.Brier)
              mean.MAE2  <- sum(RAWData.Brier$MAE2, na.rm = TRUE)/nrow(RAWData.Brier)
              mean.Brier1 <- sum(RAWData.Brier$Brier1, na.rm = TRUE)/nrow(RAWData.Brier)
              mean.Brier2 <- sum(RAWData.Brier$Brier2, na.rm = TRUE)/nrow(RAWData.Brier)
              mean.Brier[j, 1] <- mean.Brier1
              mean.Brier[j, 2] <- mean.Brier2
              mean.MAE[j, 1] <- mean.MAE1
              mean.MAE[j, 2] <- mean.MAE2
            }
            
            
          }
          
          Brier.cv[[t]] <- mean.Brier
          MAE.cv[[t]] <- mean.MAE
        } else {
          Surv <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 2))
          colnames(Surv) <- c("ID", "Surv")
          Surv$ID <- val.cdata[, ID]
          Gs <- summary(fitKM, times = landmark.time)$surv
          mean.Brier <- matrix(NA, nrow = length(horizon.time), ncol = 1)
          mean.MAE <- matrix(NA, nrow = length(horizon.time), ncol = 1)
          for (j in 1:length(horizon.time)) {
            fitKM.horizon <- try(summary(fitKM, times = horizon.time[j]), silent = TRUE)
            if ('try-error' %in% class(fitKM.horizon)) {
              mean.Brier[j, 1] <- NA
              mean.MAE[j, 1] <- NA
            } else {
              Gu <- fitKM.horizon$surv
              ## true counting process
              N1 <- vector()
              Gt <- vector()
              W.IPCW <- vector()
              for (i in 1:nrow(Surv)) {
                if (val.cdata[i, surv.var[1]] <= horizon.time[j] && val.cdata[i, surv.var[2]] == 1) {
                  N1[i] <- 1
                  Gt[i] <- summary(fitKM, times = val.cdata[i, surv.var[1]])$surv
                } else {
                  N1[i] <- 0
                  Gt[i] <- NA
                }
                
                if (val.cdata[i, surv.var[1]] > horizon.time[j]) {
                  W.IPCW[i] <- 1/(Gu/Gs)
                } else if (val.cdata[i, surv.var[1]] <= horizon.time[j] && val.cdata[i, surv.var[2]] == 1) {
                  W.IPCW[i] <- 1/(Gt[i]/Gs)
                } else {
                  W.IPCW[i] <- NA
                }
              }
              ## extract estimated Survival probability
              for (k in 1:nrow(Surv)) {
                Surv[k, 2] <- survfit$Pred[[k]][j, 2]
              }
              
              RAWData.Brier <- data.frame(Surv, N1, W.IPCW)
              colnames(RAWData.Brier)[1:2] <- c("ID", "Surv")
              RAWData.Brier$Brier <- RAWData.Brier$W.IPCW*
                abs(1 - RAWData.Brier$Surv - RAWData.Brier$N1)^2
              RAWData.Brier$MAE <- RAWData.Brier$W.IPCW*
                abs(1 - RAWData.Brier$Surv - RAWData.Brier$N1)^1
              mean.MAE[j, 1]  <- sum(RAWData.Brier$MAE, na.rm = TRUE)/nrow(RAWData.Brier)
              mean.Brier[j, 1] <- sum(RAWData.Brier$Brier, na.rm = TRUE)/nrow(RAWData.Brier)
            }
            
            
          }
          Brier.cv[[t]] <- mean.Brier
          MAE.cv[[t]] <- mean.MAE
        }
        
      }
      
    }
    writeLines(paste0("The ", t, " th validation is done!"))
    
  }
  result <- list(n.cv = n.cv, Brier.cv = Brier.cv, MAE.cv = MAE.cv, landmark.time = landmark.time,
                 horizon.time = horizon.time, method = method, quadpoint = quadpoint, 
                 CompetingRisk = CompetingRisk, seed = seed)
  class(result) <- "PEJMMLSM"
  
  return(result)
  
}