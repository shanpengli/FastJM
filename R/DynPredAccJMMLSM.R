##' @title Dynamic prediction accuracy metrics for JMMLSM models
##' @name DynPredAccJMMLSM
##' @aliases DynPredAccJMMLSM
##' @description
##' Computes dynamic prediction accuracy measures for fitted joint models of class
##' \code{JMMLSM} using grouped cross-validation.
##'
##' At each cross-validation fold, the function refits the joint model on the
##' training set, obtains subject-specific dynamic predictions on the validation
##' set at the requested horizon times, and evaluates the requested prediction
##' accuracy metrics.
##'
##' Supported metrics include:
##' \describe{
##'   \item{\code{"AUC"}}{time-dependent area under the ROC curve.}
##'   \item{\code{"Cindex"}}{concordance index for dynamic predictions.}
##'   \item{\code{"Brier"}}{inverse-probability-of-censoring weighted Brier score.}
##'   \item{\code{"MAE"}}{inverse-probability-of-censoring weighted mean absolute error.}
##'   \item{\code{"MAEQ"}}{quantile-based calibration summaries comparing empirical
##'   and predicted risks or survival probabilities.}
##' }
##'
##' For competing-risks models, event-specific cumulative incidence predictions
##' are evaluated. For single-event models, conditional survival probabilities
##' are evaluated.
##'
##' @param seed A numeric value used to set the random seed for cross-validation.
##'   Default is \code{100}.
##' @param object A fitted object of class \code{JMMLSM}.
##' @param landmark.time A numeric value specifying the landmark time at which
##'   dynamic prediction begins.
##' @param horizon.time A numeric vector of future times at which predicted
##'   probabilities are to be evaluated.
##' @param obs.time A character string specifying the longitudinal time variable
##'   in \code{object$ydata}.
##' @param method A character string specifying the approximation method used in
##'   dynamic prediction. Must be one of \code{"Laplace"} or \code{"GH"}.
##' @param quadpoint The number of quadrature points used when
##'   \code{method = "GH"}. If \code{NULL}, the value stored in \code{object} is
##'   used.
##' @param maxiter The maximum number of iterations allowed when refitting the
##'   model within each cross-validation fold. Default is \code{1000}.
##' @param n.cv The number of cross-validation folds. Default is \code{3}.
##' @param quantile.width A numeric value specifying the width of the quantile
##'   groups used for \code{"MAEQ"} summaries. Default is \code{0.25}. The
##'   reciprocal of \code{quantile.width} must be an integer.
##' @param opt A character string specifying the optimizer used when refitting
##'   the model in each fold. Default is \code{"nlminb"}.
##' @param initial.para Logical; if \code{TRUE}, parameter estimates from
##'   \code{object} are used as initial values in cross-validation refits.
##'   Default is \code{FALSE}.
##' @param LOCF Logical; if \code{TRUE}, the
##'   last-observation-carried-forward approach is applied for time-dependent
##'   survival covariates during prediction. Default is \code{FALSE}.
##' @param LOCFcovariate A character vector specifying the time-dependent
##'   survival covariates to be updated by last observation carried forward when
##'   \code{LOCF = TRUE}. Default is \code{NULL}.
##' @param clongdata A long-format data frame containing the time-dependent
##'   survival covariates used when \code{LOCF = TRUE}. Default is \code{NULL}.
##' @param metrics A character vector specifying which prediction accuracy
##'   metrics to compute. Supported choices are \code{"AUC"},
##'   \code{"Cindex"}, \code{"Brier"}, \code{"MAE"}, and \code{"MAEQ"}.
##'   Default is \code{c("AUC", "Cindex", "Brier", "MAE", "MAEQ")}.
##' @param ... Further arguments passed to or from other methods.
##'
##' @return
##' An object of class \code{DynPredAccJMMLSM}, returned as a list containing:
##' \describe{
##'   \item{\code{n.cv}}{The number of cross-validation folds.}
##'   \item{\code{landmark.time}}{The landmark time used for dynamic prediction.}
##'   \item{\code{horizon.time}}{The vector of horizon times used for evaluation.}
##'   \item{\code{method}}{The dynamic prediction approximation method used.}
##'   \item{\code{quadpoint}}{The number of quadrature points used, if applicable.}
##'   \item{\code{CompetingRisk}}{Logical; indicates whether the fitted model
##'   accounts for competing risks.}
##'   \item{\code{seed}}{The random seed used for cross-validation.}
##'   \item{\code{metrics}}{The requested evaluation metrics.}
##'   \item{\code{quantile.width}}{The width of quantile groups used for
##'   quantile-based calibration summaries.}
##'   \item{\code{AUC.cv}}{A list of fold-specific time-dependent AUC estimates,
##'   returned when \code{"AUC"} is requested. For competing-risks models, each
##'   fold contains a matrix with one column per event type.}
##'   \item{\code{Cindex.cv}}{A list of fold-specific concordance index estimates,
##'   returned when \code{"Cindex"} is requested. For competing-risks models, each
##'   fold contains a matrix with one column per event type.}
##'   \item{\code{Brier.cv}}{A list of fold-specific inverse-probability-of-
##'   censoring weighted Brier scores, returned when \code{"Brier"} is requested.
##'   For competing-risks models, each fold contains a matrix with one column per
##'   event type.}
##'   \item{\code{MAE.cv}}{A list of fold-specific inverse-probability-of-
##'   censoring weighted mean absolute errors, returned when \code{"MAE"} is
##'   requested. For competing-risks models, each fold contains a matrix with one
##'   column per event type.}
##'   \item{\code{MAEQ.cv}}{A list of fold-specific quantile-based calibration
##'   summaries, returned when \code{"MAEQ"} is requested.}
##' }
##'
##' @details
##' The function performs grouped cross-validation at the subject level. Within
##' each fold, the model is refit on the training data, and dynamic predictions
##' are computed for the validation subjects who remain under observation beyond
##' the landmark time and who have longitudinal observations observed up to that
##' time. For single-event models, prediction accuracy is based on conditional
##' survival probabilities. For competing-risks models, prediction accuracy is
##' based on event-specific cumulative incidence functions.
##' 
##' @seealso \code{\link{JMMLSM}}, \code{\link{survfitJMMLSM}}
##'
##' @export
DynPredAccJMMLSM <- function(seed = 100, object, landmark.time = NULL,
                             horizon.time = NULL, obs.time = NULL,
                             method = c("Laplace", "GH"),
                             quadpoint = NULL, maxiter = 1000,
                             n.cv = 3,
                             quantile.width = 0.25,
                             opt = "nlminb", initial.para = FALSE,
                             LOCF = FALSE, LOCFcovariate = NULL,
                             clongdata = NULL,
                             metrics = c("AUC", "Cindex", "Brier", "MAE", "MAEQ"),
                             ...) {
  
  if (!inherits(object, "JMMLSM"))
    stop("Use only with 'JMMLSM' objects.\n")
  
  if (is.null(landmark.time))
    stop("Please specify the landmark.time for dynamic prediction.")
  
  method <- match.arg(method)
  
  if (is.null(horizon.time) || !is.vector(horizon.time))
    stop("horizon.time must be vector typed.")
  
  if (is.null(obs.time)) {
    stop("Please specify a vector that represents the time variable from ydatanew.")
  } else {
    if (!obs.time %in% colnames(object$ydata)) {
      stop(paste0(obs.time, " is not found in ynewdata."))
    }
  }
  
  if (is.null(quadpoint)) {
    quadpoint <- object$quadpoint
  }
  
  allowed.metrics <- c("AUC", "Cindex", "Brier", "MAE", "MAEQ")
  if (length(metrics) < 1 || any(!metrics %in% allowed.metrics)) {
    stop("Please choose metrics from: 'AUC', 'Cindex', 'Brier', 'MAE', 'MAEQ'.")
  }
  
  if ("MAEQ" %in% metrics) {
    groups <- 1 / quantile.width
    if (floor(groups) != groups) {
      stop("The reciprocal of quantile.width must be an integer.")
    }
  } else {
    groups <- NULL
  }
  
  CompetingRisk <- object$CompetingRisk
  set.seed(seed)
  
  cdata <- object$cdata
  ydata <- object$ydata
  long.formula <- object$LongitudinalSubmodelmean
  surv.formula <- object$SurvivalSubmodel
  surv.var <- all.vars(surv.formula)
  
  variance.formula <- as.formula(
    paste("", object$LongitudinalSubmodelvariance[3], sep = "~")
  )
  
  random <- all.vars(object$random)
  ID <- random[length(random)]
  
  New.surv.formula.out <- paste0(
    "survival::Surv(", surv.var[1], ",", surv.var[2], "==0)"
  )
  New.surv.formula <- as.formula(
    paste(New.surv.formula.out, 1, sep = "~")
  )
  
  if (initial.para) {
    initial.para <- list(
      beta = object$beta,
      tau = object$tau,
      gamma1 = object$gamma1,
      gamma2 = object$gamma2,
      alpha1 = object$alpha1,
      alpha2 = object$alpha2,
      vee1 = object$vee1,
      vee2 = object$vee2,
      Sig = object$Sig
    )
  } else {
    initial.para <- NULL
  }
  
  AUC.cv <- if ("AUC" %in% metrics) list() else NULL
  Cindex.cv <- if ("Cindex" %in% metrics) list() else NULL
  Brier.cv <- if ("Brier" %in% metrics) list() else NULL
  MAE.cv <- if ("MAE" %in% metrics) list() else NULL
  MAEQ.cv <- if ("MAEQ" %in% metrics) list() else NULL
  
  folds <- caret::groupKFold(c(1:nrow(cdata)), k = n.cv)
  
  for (t in seq_len(n.cv)) {
    
    train.cdata <- cdata[folds[[t]], ]
    train.ydata <- ydata[ydata[, ID] %in% train.cdata[, ID], ]
    
    fit <- try(
      JMMLSM(
        cdata = train.cdata,
        ydata = train.ydata,
        long.formula = long.formula,
        surv.formula = surv.formula,
        variance.formula = variance.formula,
        random = object$random,
        control = JMMLSM_control(
          maxiter = maxiter,
          tol = object$tol,
          quadpoint = quadpoint,
          initial.para = initial.para,
          method = object$method,
          opt = opt
        )
      ),
      silent = TRUE
    )
    
    if ("try-error" %in% class(fit)) {
      writeLines(paste0("Error occurred in the ", t, " th training!"))
      
      if (!is.null(AUC.cv)) AUC.cv[[t]] <- NULL
      if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- NULL
      if (!is.null(Brier.cv)) Brier.cv[[t]] <- NULL
      if (!is.null(MAE.cv)) MAE.cv[[t]] <- NULL
      if (!is.null(MAEQ.cv)) MAEQ.cv[[t]] <- NULL
      
      next
    }
    
    if (fit$iter == maxiter) {
      writeLines(paste0("The ", t, " th training reached maxiter and was skipped."))
      
      if (!is.null(AUC.cv)) AUC.cv[[t]] <- NULL
      if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- NULL
      if (!is.null(Brier.cv)) Brier.cv[[t]] <- NULL
      if (!is.null(MAE.cv)) MAE.cv[[t]] <- NULL
      if (!is.null(MAEQ.cv)) MAEQ.cv[[t]] <- NULL
      
      next
    }
    
    val.cdata <- cdata[-folds[[t]], ]
    val.ydata <- ydata[ydata[, ID] %in% val.cdata[, ID], ]
    
    if ("Brier" %in% metrics || "MAE" %in% metrics) {
      fitKM <- survival::survfit(New.surv.formula, data = val.cdata)
      Gs <- summary(fitKM, times = landmark.time)$surv
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
      survfitJMMLSM(
        fit,
        ynewdata = val.ydata,
        cnewdata = val.cdata,
        u = horizon.time,
        method = method,
        Last.time = landmark.time,
        obs.time = obs.time,
        quadpoint = quadpoint,
        LOCF = LOCF,
        LOCFcovariate = LOCFcovariate,
        clongdata = val.clongdata
      ),
      silent = TRUE
    )
    
    if ("try-error" %in% class(survfit)) {
      writeLines(paste0("Error occurred in the ", t, " th validation!"))
      
      if (!is.null(AUC.cv)) AUC.cv[[t]] <- NULL
      if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- NULL
      if (!is.null(Brier.cv)) Brier.cv[[t]] <- NULL
      if (!is.null(MAE.cv)) MAE.cv[[t]] <- NULL
      if (!is.null(MAEQ.cv)) MAEQ.cv[[t]] <- NULL
      
      next
    }
    
    if (CompetingRisk) {
      
      CIF <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 3))
      colnames(CIF) <- c("ID", "CIF1", "CIF2")
      CIF$ID <- val.cdata[, ID]
      CIF$time <- val.cdata[, surv.var[1]]
      CIF$status <- val.cdata[, surv.var[2]]
      
      mean.AUC <- if ("AUC" %in% metrics) {
        matrix(NA, nrow = length(horizon.time), ncol = 2,
               dimnames = list(horizon.time, c("event1", "event2")))
      } else NULL
      
      mean.Cindex <- if ("Cindex" %in% metrics) {
        matrix(NA, nrow = length(horizon.time), ncol = 2,
               dimnames = list(horizon.time, c("event1", "event2")))
      } else NULL
      
      mean.Brier <- if ("Brier" %in% metrics) {
        matrix(NA, nrow = length(horizon.time), ncol = 2,
               dimnames = list(horizon.time, c("event1", "event2")))
      } else NULL
      
      mean.MAE <- if ("MAE" %in% metrics) {
        matrix(NA, nrow = length(horizon.time), ncol = 2,
               dimnames = list(horizon.time, c("event1", "event2")))
      } else NULL
      
      if ("MAEQ" %in% metrics) {
        AllCIF1 <- list()
        AllCIF2 <- list()
      }
      
      for (j in seq_along(horizon.time)) {
        
        for (k in seq_len(nrow(CIF))) {
          CIF[k, 2] <- survfit$Pred[[k]][j, 2]
          CIF[k, 3] <- survfit$Pred[[k]][j, 3]
        }
        
        if ("AUC" %in% metrics) {
          
          ROC1 <- timeROC(
            T = CIF$time,
            delta = CIF$status,
            weighting = "marginal",
            marker = CIF$CIF1,
            cause = 1,
            times = horizon.time[j]
          )
          
          mean.AUC[j, 1] <- ROC1$AUC_1[2]
          
          CIF$status2 <- ifelse(
            CIF$status == 2, 1,
            ifelse(CIF$status == 1, 2, 0)
          )
          
          ROC2 <- timeROC(
            T = CIF$time,
            delta = CIF$status2,
            weighting = "marginal",
            marker = CIF$CIF2,
            cause = 1,
            times = horizon.time[j]
          )
          
          mean.AUC[j, 2] <- ROC2$AUC_1[2]
        }
        
        if ("Cindex" %in% metrics) {
          mean.Cindex[j, 1] <- CindexCR(CIF$time, CIF$status, CIF$CIF1, 1)
          mean.Cindex[j, 2] <- CindexCR(CIF$time, CIF$status, CIF$CIF2, 2)
        }
        
        if ("Brier" %in% metrics || "MAE" %in% metrics) {
          
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
            
            for (i in seq_len(nrow(CIF))) {
              
              N1[i] <- as.numeric(
                val.cdata[i, surv.var[1]] <= horizon.time[j] &&
                  val.cdata[i, surv.var[2]] == 1
              )
              
              N2[i] <- as.numeric(
                val.cdata[i, surv.var[1]] <= horizon.time[j] &&
                  val.cdata[i, surv.var[2]] == 2
              )
              
              if (val.cdata[i, surv.var[1]] <= horizon.time[j] &&
                  val.cdata[i, surv.var[2]] != 0) {
                Gt[i] <- summary(fitKM, times = val.cdata[i, surv.var[1]])$surv
              } else {
                Gt[i] <- NA
              }
              
              if (val.cdata[i, surv.var[1]] > horizon.time[j]) {
                W.IPCW[i] <- 1 / (Gu / Gs)
              } else if (val.cdata[i, surv.var[1]] <= horizon.time[j] &&
                         val.cdata[i, surv.var[2]] != 0) {
                W.IPCW[i] <- 1 / (Gt[i] / Gs)
              } else {
                W.IPCW[i] <- NA
              }
            }
            
            RAWData.Brier <- data.frame(CIF, N1, N2, W.IPCW)
            
            if ("Brier" %in% metrics) {
              RAWData.Brier$Brier1 <- RAWData.Brier$W.IPCW *
                abs(RAWData.Brier$CIF1 - RAWData.Brier$N1)^2
              RAWData.Brier$Brier2 <- RAWData.Brier$W.IPCW *
                abs(RAWData.Brier$CIF2 - RAWData.Brier$N2)^2
              
              mean.Brier[j, 1] <- sum(RAWData.Brier$Brier1, na.rm = TRUE) /
                nrow(RAWData.Brier)
              mean.Brier[j, 2] <- sum(RAWData.Brier$Brier2, na.rm = TRUE) /
                nrow(RAWData.Brier)
            }
            
            if ("MAE" %in% metrics) {
              RAWData.Brier$MAE1 <- RAWData.Brier$W.IPCW *
                abs(RAWData.Brier$CIF1 - RAWData.Brier$N1)
              RAWData.Brier$MAE2 <- RAWData.Brier$W.IPCW *
                abs(RAWData.Brier$CIF2 - RAWData.Brier$N2)
              
              mean.MAE[j, 1] <- sum(RAWData.Brier$MAE1, na.rm = TRUE) /
                nrow(RAWData.Brier)
              mean.MAE[j, 2] <- sum(RAWData.Brier$MAE2, na.rm = TRUE) /
                nrow(RAWData.Brier)
            }
          }
        }
        
        if ("MAEQ" %in% metrics) {
          
          quant1 <- quantile(CIF$CIF1, probs = seq(0, 1, by = quantile.width))
          EmpiricalCIF1 <- rep(NA, groups)
          PredictedCIF1 <- rep(NA, groups)
          
          for (i in seq_len(groups)) {
            
            subquant <- CIF[
              CIF$CIF1 > quant1[i] & CIF$CIF1 <= quant1[i + 1],
              1:2
            ]
            
            quantsubdata <- val.cdata[
              val.cdata[, ID] %in% subquant$ID,
              surv.var,
              drop = FALSE
            ]
            
            quantsubCIF <- GetEmpiricalCIF(
              data = quantsubdata,
              time = surv.var[1],
              status = surv.var[2]
            )
            
            quantsubRisk1 <- quantsubCIF$H1
            
            ii <- 1
            while (ii <= nrow(quantsubRisk1)) {
              if (quantsubRisk1[ii, 1] > horizon.time[j]) {
                EmpiricalCIF1[i] <- if (ii >= 2) quantsubRisk1[ii - 1, 4] else 0
                break
              } else {
                ii <- ii + 1
              }
            }
            
            if (is.na(EmpiricalCIF1[i])) {
              EmpiricalCIF1[i] <- if (nrow(quantsubRisk1) == 0) {
                0
              } else {
                quantsubRisk1[nrow(quantsubRisk1), 4]
              }
            }
            
            PredictedCIF1[i] <- mean(subquant$CIF1)
          }
          
          AllCIF1[[j]] <- data.frame(EmpiricalCIF1, PredictedCIF1)
          
          quant2 <- quantile(CIF$CIF2, probs = seq(0, 1, by = quantile.width))
          EmpiricalCIF2 <- rep(NA, groups)
          PredictedCIF2 <- rep(NA, groups)
          
          for (i in seq_len(groups)) {
            
            subquant <- CIF[
              CIF$CIF2 > quant2[i] & CIF$CIF2 <= quant2[i + 1],
              c(1, 3)
            ]
            
            quantsubdata <- val.cdata[
              val.cdata[, ID] %in% subquant$ID,
              surv.var,
              drop = FALSE
            ]
            
            quantsubCIF <- GetEmpiricalCIF(
              data = quantsubdata,
              time = surv.var[1],
              status = surv.var[2]
            )
            
            quantsubRisk2 <- quantsubCIF$H2
            
            ii <- 1
            while (ii <= nrow(quantsubRisk2)) {
              if (quantsubRisk2[ii, 1] > horizon.time[j]) {
                EmpiricalCIF2[i] <- if (ii >= 2) quantsubRisk2[ii - 1, 4] else 0
                break
              } else {
                ii <- ii + 1
              }
            }
            
            if (is.na(EmpiricalCIF2[i])) {
              EmpiricalCIF2[i] <- if (nrow(quantsubRisk2) == 0) {
                0
              } else {
                quantsubRisk2[nrow(quantsubRisk2), 4]
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
      
      mean.AUC <- if ("AUC" %in% metrics) {
        matrix(NA, nrow = length(horizon.time), ncol = 1,
               dimnames = list(horizon.time, "event"))
      } else NULL
      
      mean.Cindex <- if ("Cindex" %in% metrics) {
        matrix(NA, nrow = length(horizon.time), ncol = 1,
               dimnames = list(horizon.time, "event"))
      } else NULL
      
      mean.Brier <- if ("Brier" %in% metrics) {
        matrix(NA, nrow = length(horizon.time), ncol = 1,
               dimnames = list(horizon.time, "event"))
      } else NULL
      
      mean.MAE <- if ("MAE" %in% metrics) {
        matrix(NA, nrow = length(horizon.time), ncol = 1,
               dimnames = list(horizon.time, "event"))
      } else NULL
      
      if ("MAEQ" %in% metrics) {
        AllSurv <- list()
      }
      
      for (j in seq_along(horizon.time)) {
        
        for (k in seq_len(nrow(Surv))) {
          Surv[k, 2] <- survfit$Pred[[k]][j, 2]
        }
        
        if ("AUC" %in% metrics) {
          
          ROC <- timeROC(
            T = Surv$time,
            delta = Surv$status,
            weighting = "marginal",
            marker = -Surv$Surv,
            cause = 1,
            times = horizon.time[j]
          )
          
          mean.AUC[j, 1] <- ROC$AUC[2]
        }
        
        if ("Cindex" %in% metrics) {
          mean.Cindex[j, 1] <- CindexCR(
            Surv$time,
            Surv$status,
            1 - Surv$Surv,
            1
          )
        }
        
        if ("Brier" %in% metrics || "MAE" %in% metrics) {
          
          fitKM.horizon <- try(summary(fitKM, times = horizon.time[j]), silent = TRUE)
          
          if ("try-error" %in% class(fitKM.horizon)) {
            
            if ("Brier" %in% metrics) mean.Brier[j, 1] <- NA
            if ("MAE" %in% metrics) mean.MAE[j, 1] <- NA
            
          } else {
            
            Gu <- fitKM.horizon$surv
            
            N1 <- vector()
            Gt <- vector()
            W.IPCW <- vector()
            
            for (i in seq_len(nrow(Surv))) {
              
              if (val.cdata[i, surv.var[1]] <= horizon.time[j] &&
                  val.cdata[i, surv.var[2]] == 1) {
                N1[i] <- 1
                Gt[i] <- summary(fitKM, times = val.cdata[i, surv.var[1]])$surv
              } else {
                N1[i] <- 0
                Gt[i] <- NA
              }
              
              if (val.cdata[i, surv.var[1]] > horizon.time[j]) {
                W.IPCW[i] <- 1 / (Gu / Gs)
              } else if (val.cdata[i, surv.var[1]] <= horizon.time[j] &&
                         val.cdata[i, surv.var[2]] == 1) {
                W.IPCW[i] <- 1 / (Gt[i] / Gs)
              } else {
                W.IPCW[i] <- NA
              }
            }
            
            RAWData.Brier <- data.frame(Surv, N1, W.IPCW)
            
            if ("Brier" %in% metrics) {
              RAWData.Brier$Brier <- RAWData.Brier$W.IPCW *
                abs(1 - RAWData.Brier$Surv - RAWData.Brier$N1)^2
              
              mean.Brier[j, 1] <- sum(RAWData.Brier$Brier, na.rm = TRUE) /
                nrow(RAWData.Brier)
            }
            
            if ("MAE" %in% metrics) {
              RAWData.Brier$MAE <- RAWData.Brier$W.IPCW *
                abs(1 - RAWData.Brier$Surv - RAWData.Brier$N1)
              
              mean.MAE[j, 1] <- sum(RAWData.Brier$MAE, na.rm = TRUE) /
                nrow(RAWData.Brier)
            }
          }
        }
        
        if ("MAEQ" %in% metrics) {
          
          quant <- quantile(Surv$Surv, probs = seq(0, 1, by = quantile.width))
          EmpiricalSurv <- rep(NA, groups)
          PredictedSurv <- rep(NA, groups)
          
          for (i in seq_len(groups)) {
            
            subquant <- Surv[
              Surv$Surv > quant[i] & Surv$Surv <= quant[i + 1],
              c(1, 2)
            ]
            
            quantsubdata <- val.cdata[
              val.cdata[, ID] %in% subquant$ID,
              surv.var,
              drop = FALSE
            ]
            
            colnames(quantsubdata) <- c("time", "status")
            
            fitKM.sub <- survival::survfit(
              survival::Surv(time, status) ~ 1,
              data = quantsubdata
            )
            
            fitKM.horizon.sub <- try(
              summary(fitKM.sub, times = horizon.time[j]),
              silent = TRUE
            )
            
            if ("try-error" %in% class(fitKM.horizon.sub)) {
              EmpiricalSurv[i] <- summary(
                fitKM.sub,
                times = max(quantsubdata$time)
              )$surv
            } else {
              tempSurv <- summary(fitKM.sub, times = horizon.time[j])$surv
              EmpiricalSurv[i] <- if (is.null(tempSurv)) 0 else tempSurv
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
  
  class(result) <- "DynPredAccJMMLSM"
  return(result)
}