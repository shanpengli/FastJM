##' @title Calculating evaluation metrics for joint models
##' @name DynPredAccjmcs
##' @aliases DynPredAccjmcs
##' @param seed a numeric value of seed to be specified for cross validation.
##' @param object object of class 'jmcs'.
##' @param landmark.time a numeric value of time for which dynamic prediction starts.
##' @param horizon.time a numeric vector of future times for which predicted probabilities are to be computed.
##' @param obs.time a character string of specifying a longitudinal time variable.
##' @param method estimation method for predicted probabilities. If \code{Laplace}, then the empirical empirical
##' estimates of random effects is used. If \code{GH}, then the pseudo-adaptive Gauss-Hermite quadrature is used.
##' @param quadpoint the number of pseudo-adaptive Gauss-Hermite quadrature points if \code{method = "GH"}.
##' @param maxiter the maximum number of iterations of the EM algorithm that the 
##' function will perform. Default is 10000.
##' @param n.cv number of folds for cross validation. Default is 3.
##' @param survinitial Fit a Cox model to obtain initial values of the parameter estimates. Default is TRUE.
##' @param quantile.width a numeric value of width of quantile to be specified. Default is 0.25.
##' @param initial.para Initial guess of parameters for cross validation. Default is FALSE.
##' @param LOCF a logical value to indicate whether the last-observation-carried-forward approach applies to prediction. 
##' If \code{TRUE}, then \code{LOCFcovariate} and \code{clongdata} must be specified to indicate 
##' which time-dependent survival covariates are included for dynamic prediction. Default is FALSE.
##' @param LOCFcovariate a vector of string with time-dependent survival covariates if \code{LOCF = TRUE}. Default is NULL.
##' @param clongdata a long format data frame where time-dependent survival covariates are incorporated. Default is NULL.
##' @param metrics a list to indicate which metric is used. 
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @seealso \code{\link{jmcs}, \link{survfitjmcs}}
##' @export
##' 
DynPredAccjmcs <- function(seed = 100, object, landmark.time = NULL, horizon.time = NULL,
                     obs.time = NULL, method = c("Laplace", "GH"),
                     quadpoint = NULL, maxiter = NULL, n.cv = 3,
                     survinitial = TRUE,
                     quantile.width = 0.25,
                     initial.para = FALSE,
                     LOCF = FALSE, LOCFcovariate = NULL, clongdata = NULL,
                     metrics = c("AUC", "Cindex", "Brier", "MAE", "MAEQ"), ...) {
  
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' xs.\n")
  
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
  
  if (is.null(maxiter)) {
    maxiter <- 10000
  }
  
  allowed.metrics <- c("AUC", "Cindex", "Brier", "MAE", "MAEQ")
  if (length(metrics) < 1 || any(!metrics %in% allowed.metrics)) {
    stop("Please choose metrics from: 'AUC', 'Cindex', 'Brier', 'MAE', 'MAEQ'.")
  }
  
  if ("MAEQ" %in% metrics) {
    groups <- 1 / quantile.width
    if (floor(groups) != groups)
      stop("The reciprocal of quantile.width must be an integer.")
  } else {
    groups <- NULL
  }
  
  CompetingRisk <- object$CompetingRisk
  set.seed(seed)
  cdata <- object$cdata
  ydata <- object$ydata
  long.formula <- object$LongitudinalSubmodel
  surv.formula <- object$SurvivalSubmodel
  surv.var <- all.vars(surv.formula)
  
  New.surv.formula.out <- paste0("survival::Surv(", surv.var[1], ",", surv.var[2], "==0)")
  New.surv.formula <- as.formula(paste(New.surv.formula.out, 1, sep = "~")) # for prediction error (PE)
  
  random <- all.vars(object$random)
  ID <- random[length(random)]
  
  if (initial.para) {
    initial.para <- list(beta = object$beta,
                         sigma = object$sigma,
                         gamma1 = object$gamma1,
                         gamma2 = object$gamma2,
                         alpha1 = object$nu1,
                         alpha2 = object$nu2,
                         Sig = object$Sig)
  } else {
    initial.para <- NULL
  }
  
  need.discrimination <- any(metrics %in% c("AUC", "Cindex"))
  need.error <- any(metrics %in% c("Brier", "MAE"))
  need.maeq <- "MAEQ" %in% metrics
  
  AUC.cv <- if ("AUC" %in% metrics) list() else NULL
  Cindex.cv <- if ("Cindex" %in% metrics) list() else NULL
  Brier.cv <- if ("Brier" %in% metrics) list() else NULL
  MAE.cv <- if ("MAE" %in% metrics) list() else NULL
  MAEQ.cv <- if ("MAEQ" %in% metrics) list() else NULL
  
  folds <- caret::groupKFold(c(1:nrow(cdata)), k = n.cv)
  
  for (t in 1:n.cv) {
    
    train.cdata <- cdata[folds[[t]], ]
    train.ydata <- ydata[ydata[, ID] %in% train.cdata[, ID], ]
    
    fit <- try(
      jmcs(cdata = train.cdata,
           ydata = train.ydata,
           long.formula = long.formula,
           surv.formula = surv.formula,
           quadpoint = quadpoint,
           random = object$random,
           survinitial = survinitial,
           opt = object$opt,
           initial.para = initial.para),
      silent = TRUE
    )
    
    if ("try-error" %in% class(fit)) {
      writeLines(paste0("Error occured in the ", t, " th training!"))
      
      if (!is.null(AUC.cv)) AUC.cv[[t]] <- NULL
      if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- NULL
      if (!is.null(Brier.cv)) Brier.cv[[t]] <- NULL
      if (!is.null(MAE.cv)) MAE.cv[[t]] <- NULL
      if (!is.null(MAEQ.cv)) MAEQ.cv[[t]] <- NULL
      
    } else if (fit$iter == maxiter) {
      
      if (!is.null(AUC.cv)) AUC.cv[[t]] <- NULL
      if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- NULL
      if (!is.null(Brier.cv)) Brier.cv[[t]] <- NULL
      if (!is.null(MAE.cv)) MAE.cv[[t]] <- NULL
      if (!is.null(MAEQ.cv)) MAEQ.cv[[t]] <- NULL
      
    } else {
      
      val.cdata <- cdata[-folds[[t]], ]
      val.ydata <- ydata[ydata[, ID] %in% val.cdata[, ID], ]
      
      if (need.error) {
        ## fit a Kaplan-Meier estimator on original validation cdata
        fitKM <- survival::survfit(New.surv.formula, data = val.cdata)
      } else {
        fitKM <- NULL
      }
      
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
      
      survfit <- try(
        survfitjmcs(fit,
                    ynewdata = val.ydata,
                    cnewdata = val.cdata,
                    u = horizon.time,
                    method = method,
                    Last.time = landmark.time,
                    obs.time = obs.time,
                    quadpoint = quadpoint,
                    LOCF = LOCF,
                    LOCFcovariate = LOCFcovariate,
                    clongdata = val.clongdata),
        silent = TRUE
      )
      
      if ("try-error" %in% class(survfit)) {
        writeLines(paste0("Error occured in the ", t, " th validation!"))
        
        if (!is.null(AUC.cv)) AUC.cv[[t]] <- NULL
        if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- NULL
        if (!is.null(Brier.cv)) Brier.cv[[t]] <- NULL
        if (!is.null(MAE.cv)) MAE.cv[[t]] <- NULL
        if (!is.null(MAEQ.cv)) MAEQ.cv[[t]] <- NULL
        
      } else {
        
        if (CompetingRisk) {
          
          CIF <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 3))
          colnames(CIF) <- c("ID", "CIF1", "CIF2")
          CIF$ID <- val.cdata[, ID]
          CIF$time <- val.cdata[, surv.var[1]]
          CIF$status <- val.cdata[, surv.var[2]]
          
          mean.AUC <- if ("AUC" %in% metrics) matrix(NA, nrow = length(horizon.time), ncol = 2) else NULL
          mean.Cindex <- if ("Cindex" %in% metrics) matrix(NA, nrow = length(horizon.time), ncol = 2) else NULL
          mean.Brier <- if ("Brier" %in% metrics) matrix(NA, nrow = length(horizon.time), ncol = 2) else NULL
          mean.MAE <- if ("MAE" %in% metrics) matrix(NA, nrow = length(horizon.time), ncol = 2) else NULL
          
          if (need.maeq) {
            AllCIF1 <- list()
            AllCIF2 <- list()
          } else {
            AllCIF1 <- NULL
            AllCIF2 <- NULL
          }
          
          if (need.error) {
            Gs <- summary(fitKM, times = landmark.time)$surv # for prediction error (PE)
          }
          
          for (j in 1:length(horizon.time)) {
            
            ## extract estimated CIF
            for (k in 1:nrow(CIF)) {
              CIF[k, 2] <- survfit$Pred[[k]][j, 2]
              CIF[k, 3] <- survfit$Pred[[k]][j, 3]
            }
            
            if ("AUC" %in% metrics) {
              ROC <- timeROC(T = CIF$time, delta = CIF$status,
                                      weighting = "marginal",
                                      marker = CIF$CIF1, cause = 1,
                                      times = horizon.time[j])
              
              mean.AUC[j, 1] <- ROC$AUC_1[2]
              
              CIF$status2 <- ifelse(CIF$status == 2, 1,
                                    ifelse(CIF$status == 1, 2, 0))
              
              ROC <- timeROC(T = CIF$time, delta = CIF$status2,
                                      weighting = "marginal",
                                      marker = CIF$CIF2, cause = 1,
                                      times = horizon.time[j])
              
              mean.AUC[j, 2] <- ROC$AUC_1[2]
            }
            
            if ("Cindex" %in% metrics) {
              mean.Cindex[j, 1] <- CindexCR(CIF$time, CIF$status, CIF$CIF1, 1)
              mean.Cindex[j, 2] <- CindexCR(CIF$time, CIF$status, CIF$CIF2, 2)
            }
            
            if (need.error) {
              fitKM.horizon <- try(summary(fitKM, times = horizon.time[j]), silent = TRUE)
              
              if ("try-error" %in% class(fitKM.horizon)) {
                if ("Brier" %in% metrics) {
                  mean.Brier[j, 1] <- NA
                  mean.Brier[j, 2] <- NA
                }
                if ("MAE" %in% metrics) {
                  mean.MAE[j, 1] <- NA
                  mean.MAE[j, 2] <- NA
                }
              } else {
                Gu <- fitKM.horizon$surv
                
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
                
                RAWData.Brier <- data.frame(CIF, N1, N2, W.IPCW)
                colnames(RAWData.Brier)[1:3] <- c("ID", "CIF1", "CIF2")
                
                if ("Brier" %in% metrics) {
                  RAWData.Brier$Brier1 <- RAWData.Brier$W.IPCW *
                    abs(RAWData.Brier$CIF1 - RAWData.Brier$N1)^2
                  RAWData.Brier$Brier2 <- RAWData.Brier$W.IPCW *
                    abs(RAWData.Brier$CIF2 - RAWData.Brier$N2)^2
                  
                  mean.Brier[j, 1] <- sum(RAWData.Brier$Brier1, na.rm = TRUE) / nrow(RAWData.Brier)
                  mean.Brier[j, 2] <- sum(RAWData.Brier$Brier2, na.rm = TRUE) / nrow(RAWData.Brier)
                }
                
                if ("MAE" %in% metrics) {
                  RAWData.Brier$MAE1 <- RAWData.Brier$W.IPCW *
                    abs(RAWData.Brier$CIF1 - RAWData.Brier$N1)^1
                  RAWData.Brier$MAE2 <- RAWData.Brier$W.IPCW *
                    abs(RAWData.Brier$CIF2 - RAWData.Brier$N2)^1
                  
                  mean.MAE[j, 1] <- sum(RAWData.Brier$MAE1, na.rm = TRUE) / nrow(RAWData.Brier)
                  mean.MAE[j, 2] <- sum(RAWData.Brier$MAE2, na.rm = TRUE) / nrow(RAWData.Brier)
                }
              }
            }
            
            if (need.maeq) {
              quant1 <- quantile(CIF$CIF1, probs = seq(0, 1, by = quantile.width))
              EmpiricalCIF1 <- rep(NA, groups)
              PredictedCIF1 <- rep(NA, groups)
              
              for (i in 1:groups) {
                subquant <- CIF[CIF$CIF1 > quant1[i] &
                                  CIF$CIF1 <= quant1[i + 1], 1:2]
                quantsubdata <- val.cdata[val.cdata[, ID] %in% subquant$ID, surv.var, drop = FALSE]
                
                quantsubCIF <- GetEmpiricalCIF(data = quantsubdata,
                                               time = surv.var[1],
                                               status = surv.var[2])
                
                quantsubRisk1 <- quantsubCIF$H1
                ii <- 1
                while (ii <= nrow(quantsubRisk1)) {
                  if (quantsubRisk1[ii, 1] > horizon.time[j]) {
                    if (ii >= 2) {
                      EmpiricalCIF1[i] <- quantsubRisk1[ii - 1, 4]
                    } else {
                      EmpiricalCIF1[i] <- 0
                    }
                    break
                  } else {
                    ii <- ii + 1
                  }
                }
                if (is.na(EmpiricalCIF1[i])) {
                  if (nrow(quantsubRisk1) == 0) {
                    EmpiricalCIF1[i] <- 0
                  } else {
                    EmpiricalCIF1[i] <- quantsubRisk1[nrow(quantsubRisk1), 4]
                  }
                }
                PredictedCIF1[i] <- mean(subquant$CIF1)
              }
              
              AllCIF1[[j]] <- data.frame(EmpiricalCIF1, PredictedCIF1)
              
              quant2 <- quantile(CIF$CIF2, probs = seq(0, 1, by = quantile.width))
              EmpiricalCIF2 <- rep(NA, groups)
              PredictedCIF2 <- rep(NA, groups)
              
              for (i in 1:groups) {
                subquant <- CIF[CIF$CIF2 > quant2[i] &
                                  CIF$CIF2 <= quant2[i + 1], c(1, 3)]
                quantsubdata <- val.cdata[val.cdata[, ID] %in% subquant$ID, surv.var, drop = FALSE]
                
                quantsubCIF <- GetEmpiricalCIF(data = quantsubdata,
                                               time = surv.var[1],
                                               status = surv.var[2])
                
                quantsubRisk2 <- quantsubCIF$H2
                ii <- 1
                while (ii <= nrow(quantsubRisk2)) {
                  if (quantsubRisk2[ii, 1] > horizon.time[j]) {
                    if (ii >= 2) {
                      EmpiricalCIF2[i] <- quantsubRisk2[ii - 1, 4]
                    } else {
                      EmpiricalCIF2[i] <- 0
                    }
                    break
                  } else {
                    ii <- ii + 1
                  }
                }
                if (is.na(EmpiricalCIF2[i])) {
                  if (nrow(quantsubRisk2) == 0) {
                    EmpiricalCIF2[i] <- 0
                  } else {
                    EmpiricalCIF2[i] <- quantsubRisk2[nrow(quantsubRisk2), 4]
                  }
                }
                PredictedCIF2[i] <- mean(subquant$CIF2)
              }
              
              AllCIF2[[j]] <- data.frame(EmpiricalCIF2, PredictedCIF2)
            }
          }
          
          if (!is.null(AUC.cv)) AUC.cv[[t]] <- mean.AUC
          if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- mean.Cindex
          if (!is.null(Brier.cv)) Brier.cv[[t]] <- mean.Brier
          if (!is.null(MAE.cv)) MAE.cv[[t]] <- mean.MAE
          
          if (!is.null(MAEQ.cv)) {
            names(AllCIF1) <- names(AllCIF2) <- horizon.time
            MAEQ.cv[[t]] <- list(AllCIF1 = AllCIF1, AllCIF2 = AllCIF2)
          }
          
        } else {
          
          Surv <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 2))
          colnames(Surv) <- c("ID", "Surv")
          Surv$ID <- val.cdata[, ID]
          Surv$time <- val.cdata[, surv.var[1]]
          Surv$status <- val.cdata[, surv.var[2]]
          
          mean.AUC <- if ("AUC" %in% metrics) matrix(NA, nrow = length(horizon.time), ncol = 1) else NULL
          mean.Cindex <- if ("Cindex" %in% metrics) matrix(NA, nrow = length(horizon.time), ncol = 1) else NULL
          mean.Brier <- if ("Brier" %in% metrics) matrix(NA, nrow = length(horizon.time), ncol = 1) else NULL
          mean.MAE <- if ("MAE" %in% metrics) matrix(NA, nrow = length(horizon.time), ncol = 1) else NULL
          
          if (need.maeq) {
            AllSurv <- list()
          } else {
            AllSurv <- NULL
          }
          
          if (need.error) {
            Gs <- summary(fitKM, times = landmark.time)$surv
          }
          
          for (j in 1:length(horizon.time)) {
            
            ## extract estimated Survival probability
            for (k in 1:nrow(Surv)) {
              Surv[k, 2] <- survfit$Pred[[k]][j, 2]
            }
            
            if ("AUC" %in% metrics) {
              ROC <- timeROC(T = Surv$time, delta = Surv$status,
                                      weighting = "marginal",
                                      marker = -Surv$Surv, cause = 1,
                                      times = horizon.time[j])
              
              mean.AUC[j, 1] <- ROC$AUC[2]
            }
            
            if ("Cindex" %in% metrics) {
              mean.Cindex[j, 1] <- CindexCR(Surv$time, Surv$status, -Surv$Surv, 1)
            }
            
            if (need.error) {
              fitKM.horizon <- try(summary(fitKM, times = horizon.time[j]), silent = TRUE)
              
              if ("try-error" %in% class(fitKM.horizon)) {
                if ("Brier" %in% metrics) mean.Brier[j, 1] <- NA
                if ("MAE" %in% metrics) mean.MAE[j, 1] <- NA
              } else {
                Gu <- fitKM.horizon$surv
                
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
                
                RAWData.Brier <- data.frame(Surv, N1, W.IPCW)
                colnames(RAWData.Brier)[1:2] <- c("ID", "Surv")
                
                if ("Brier" %in% metrics) {
                  RAWData.Brier$Brier <- RAWData.Brier$W.IPCW *
                    abs(1 - RAWData.Brier$Surv - RAWData.Brier$N1)^2
                  mean.Brier[j, 1] <- sum(RAWData.Brier$Brier, na.rm = TRUE) / nrow(RAWData.Brier)
                }
                
                if ("MAE" %in% metrics) {
                  RAWData.Brier$MAE <- RAWData.Brier$W.IPCW *
                    abs(1 - RAWData.Brier$Surv - RAWData.Brier$N1)^1
                  mean.MAE[j, 1] <- sum(RAWData.Brier$MAE, na.rm = TRUE) / nrow(RAWData.Brier)
                }
              }
            }
            
            if (need.maeq) {
              quant <- quantile(Surv$Surv, probs = seq(0, 1, by = quantile.width))
              EmpiricalSurv <- rep(NA, groups)
              PredictedSurv <- rep(NA, groups)
              
              for (i in 1:groups) {
                subquant <- Surv[Surv$Surv > quant[i] &
                                   Surv$Surv <= quant[i + 1], c(1, 2)]
                quantsubdata <- val.cdata[val.cdata[, ID] %in% subquant$ID, surv.var, drop = FALSE]
                colnames(quantsubdata) <- c("time", "status")
                
                fitKM.sub <- survival::survfit(survival::Surv(time, status) ~ 1, data = quantsubdata)
                fitKM.horizon.sub <- try(summary(fitKM.sub, times = horizon.time[j]), silent = TRUE)
                
                if ("try-error" %in% class(fitKM.horizon.sub)) {
                  EmpiricalSurv[i] <- summary(fitKM.sub, times = max(quantsubdata$time))$surv
                } else {
                  tempSurv <- summary(fitKM.sub, times = horizon.time[j])$surv
                  if (is.null(tempSurv)) {
                    EmpiricalSurv[i] <- 0
                  } else {
                    EmpiricalSurv[i] <- tempSurv
                  }
                  
                }
                
                PredictedSurv[i] <- mean(subquant$Surv)
              }
              
              AllSurv[[j]] <- data.frame(EmpiricalSurv, PredictedSurv)
            }
          }
          
          if (!is.null(AUC.cv)) AUC.cv[[t]] <- mean.AUC
          if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- mean.Cindex
          if (!is.null(Brier.cv)) Brier.cv[[t]] <- mean.Brier
          if (!is.null(MAE.cv)) MAE.cv[[t]] <- mean.MAE
          
          if (!is.null(MAEQ.cv)) {
            names(AllSurv) <- horizon.time
            MAEQ.cv[[t]] <- list(AllSurv = AllSurv)
          }
        }
      }
    }
    
    writeLines(paste0("The ", t, "-th validation is done!"))
  }
  
  result <- list(
    n.cv = n.cv,
    landmark.time = landmark.time,
    horizon.time = horizon.time,
    method = method,
    quadpoint = quadpoint,
    CompetingRisk = CompetingRisk,
    seed = seed,
    metrics = metrics,
    quantile.width = quantile.width,
    AUC.cv = AUC.cv,
    Cindex.cv = Cindex.cv,
    Brier.cv = Brier.cv,
    MAE.cv = MAE.cv,
    MAEQ.cv = MAEQ.cv
  )
  
  class(result) <- "DynPredAccjmcs"
  return(result)
}
