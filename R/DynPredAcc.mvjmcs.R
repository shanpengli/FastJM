DynPredAcc.mvjmcs <- function(seed = 100,
                             object,
                             landmark.time = NULL,
                             horizon.time = NULL,
                             obs.time = NULL,
                             maxiter = 1000,
                             n.cv = 3,
                             quantile.width = 0.25,
                             opt = c("nlminb", "optim"),
                             LOCF = FALSE, LOCFcovariate = NULL,
                             clongdata = NULL,
                             metrics = c("AUC", "Cindex", 
                                         "Brier Score", "MAEQ", "MAE"),
                             cpu.cores = 1, ...) {
  
  if (!inherits(object, "mvjmcs")) {
    stop("Use only with 'mvjmcs' objects.\n")
  }
  
  if (is.null(landmark.time)) {
    stop("Please specify the landmark.time for dynamic prediction.")
  }
  
  if (is.null(horizon.time) || !is.vector(horizon.time)) {
    stop("horizon.time must be vector typed.")
  }
  
  if (is.null(obs.time)) {
    stop("Please specify a vector that represents the time variable from ydata.")
  } else {
    if (!obs.time %in% colnames(object$ydata)) {
      stop(paste0(obs.time, " is not found in object$ydata."))
    }
  }
  
  metrics[metrics == "Brier Score"] <- "Brier Score"
  
  allowed.metrics <- c("AUC", "Cindex", "Brier Score",
                       "MAE", "MAEQ")
  
  if (length(metrics) < 1 || any(!metrics %in% allowed.metrics)) {
    stop(paste0("Please choose metrics from: ",
                "'AUC', 'Cindex', 'Brier Score', 'MAE', 'MAEQ'."))
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
  
  long.formula <- object$LongitudinalSubmodel
  surv.formula <- object$SurvivalSubmodel
  surv.var <- all.vars(surv.formula)
  
  random <- object$random
  
  if (is.list(random)) {
    random.vars <- all.vars(random[[1]])
  } else {
    random.vars <- all.vars(random)
  }
  
  ID <- random.vars[length(random.vars)]
  
  if (length(ID) == 0 || is.na(ID) ||
      !ID %in% colnames(cdata) ||
      !ID %in% colnames(ydata)) {
    stop(
      paste0(
        "Could not identify the subject ID variable from object$random. ",
        "Please check object$random and make sure the ID variable appears ",
        "in both object$cdata and object$ydata."
      )
    )
  }
  
  CompetingRisk <- object$CompetingRisk
  
  New.surv.formula.out <- paste0(
    "survival::Surv(", surv.var[1], ",", surv.var[2], "==0)"
  )
  
  New.surv.formula <- as.formula(
    paste(New.surv.formula.out, 1, sep = "~")
  )
  
  AUC.cv <- if ("AUC" %in% metrics) list() else NULL
  Cindex.cv <- if ("Cindex" %in% metrics) list() else NULL
  Brier.cv <- if ("Brier Score" %in% metrics) list() else NULL
  MAE.cv <- if ("MAE" %in% metrics) list() else NULL
  MAEQ.cv <- if ("MAEQ" %in% metrics) list() else NULL
  
  failed.folds <- list()
  
  folds <- caret::groupKFold(c(1:nrow(cdata)), k = n.cv)
  
  for (t in seq_len(n.cv)) {
    
    train.cdata <- cdata[folds[[t]], , drop = FALSE]
    train.ydata <- ydata[ydata[[ID]] %in% train.cdata[[ID]], , drop = FALSE]
    
    val.cdata <- cdata[-folds[[t]], , drop = FALSE]
    val.ydata <- ydata[ydata[[ID]] %in% val.cdata[[ID]], , drop = FALSE]
    
    train.cdata <- train.cdata[order(train.cdata[[ID]]), , drop = FALSE]
    train.ydata <- train.ydata[order(train.ydata[[ID]]), , drop = FALSE]
    val.cdata <- val.cdata[order(val.cdata[[ID]]), , drop = FALSE]
    val.ydata <- val.ydata[order(val.ydata[[ID]]), , drop = FALSE]
    
    if (nrow(train.cdata) == 0 || nrow(train.ydata) == 0) {
      msg <- paste0("The ", t, " th training fold has no usable subjects.")
      writeLines(msg)
      failed.folds[[as.character(t)]] <- msg
      
      if (!is.null(AUC.cv)) AUC.cv[[t]] <- NULL
      if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- NULL
      if (!is.null(Brier.cv)) Brier.cv[[t]] <- NULL
      if (!is.null(MAE.cv)) MAE.cv[[t]] <- NULL
      if (!is.null(MAEQ.cv)) MAEQ.cv[[t]] <- NULL
      
      next
    }
    
    fit <- try(
      mvjmcs(
        cdata = train.cdata,
        ydata = train.ydata,
        long.formula = long.formula,
        surv.formula = surv.formula,
        random = random,
        control = mvjmcs_control(
          cpu.cores = cpu.cores,
          maxiter = maxiter,
          tol = object$tol,
          opt = opt
        )
      ), 
      silent = TRUE
    )
    
    if ("try-error" %in% class(fit)) {
      msg <- paste0("Error occurred in the ", t, " th training!")
      writeLines(msg)
      failed.folds[[as.character(t)]] <- msg
      
      if (!is.null(AUC.cv)) AUC.cv[[t]] <- NULL
      if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- NULL
      if (!is.null(Brier.cv)) Brier.cv[[t]] <- NULL
      if (!is.null(MAE.cv)) MAE.cv[[t]] <- NULL
      if (!is.null(MAEQ.cv)) MAEQ.cv[[t]] <- NULL
      
      next
    }
    
    if (!is.null(fit$iter) && fit$iter == maxiter) {
      msg <- paste0("The ", t, " th training reached maxiter and was skipped.")
      writeLines(msg)
      failed.folds[[as.character(t)]] <- msg
      
      if (!is.null(AUC.cv)) AUC.cv[[t]] <- NULL
      if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- NULL
      if (!is.null(Brier.cv)) Brier.cv[[t]] <- NULL
      if (!is.null(MAE.cv)) MAE.cv[[t]] <- NULL
      if (!is.null(MAEQ.cv)) MAEQ.cv[[t]] <- NULL
      
      next
    }
    
    val.cdata <- val.cdata[val.cdata[, surv.var[1]] > landmark.time,
                           ,
                           drop = FALSE]
    
    val.ydata <- val.ydata[val.ydata[[ID]] %in% val.cdata[[ID]], ,
                           drop = FALSE]
    
    val.ydata <- val.ydata[val.ydata[, obs.time] <= landmark.time,
                           ,
                           drop = FALSE
                           ]
    
    NewyID <- unique(val.ydata[[ID]])
    
    val.cdata <- val.cdata[
      val.cdata[[ID]] %in% NewyID,
      ,
      drop = FALSE
    ]
    
    val.cdata <- val.cdata[order(val.cdata[[ID]]), , drop = FALSE]
    val.ydata <- val.ydata[order(val.ydata[[ID]]), , drop = FALSE]
    
    if (nrow(val.cdata) == 0 || nrow(val.ydata) == 0) {
      msg <- paste0(
        "The ", t,
        " th validation fold has no eligible subjects after landmark filtering."
      )
      writeLines(msg)
      failed.folds[[as.character(t)]] <- msg
      
      if (!is.null(AUC.cv)) AUC.cv[[t]] <- NULL
      if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- NULL
      if (!is.null(Brier.cv)) Brier.cv[[t]] <- NULL
      if (!is.null(MAE.cv)) MAE.cv[[t]] <- NULL
      if (!is.null(MAEQ.cv)) MAEQ.cv[[t]] <- NULL
      
      next
    }
    
    if (LOCF) {
      val.clongdata <- clongdata[clongdata[, ID] %in% val.cdata[, ID], ]
    } else {
      val.clongdata <- NULL
    }
    
    survfit <- try(
      survfitmvjmcs(
        fit,
        seed = seed,
        ynewdata = val.ydata,
        cnewdata = val.cdata,
        Last.time = landmark.time,
        u = horizon.time,
        obs.time = obs.time,
        LOCF = LOCF,
        LOCFcovariate = LOCFcovariate,
        clongdata = val.clongdata
      ),
      silent = TRUE
    )
    
    if ("try-error" %in% class(survfit)) {
      msg <- paste0("Error occurred in the ", t, " th validation!")
      writeLines(msg)
      failed.folds[[as.character(t)]] <- msg
      
      if (!is.null(AUC.cv)) AUC.cv[[t]] <- NULL
      if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- NULL
      if (!is.null(Brier.cv)) Brier.cv[[t]] <- NULL
      if (!is.null(MAE.cv)) MAE.cv[[t]] <- NULL
      if (!is.null(MAEQ.cv)) MAEQ.cv[[t]] <- NULL
      
      next
    }
    
    Last.time <- survfit$Last.time
    
    if (is.null(Last.time) || nrow(Last.time) == 0) {
      msg <- paste0(
        "The ", t,
        " th validation fold produced no Last.time output."
      )
      writeLines(msg)
      failed.folds[[as.character(t)]] <- msg
      
      if (!is.null(AUC.cv)) AUC.cv[[t]] <- NULL
      if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- NULL
      if (!is.null(Brier.cv)) Brier.cv[[t]] <- NULL
      if (!is.null(MAE.cv)) MAE.cv[[t]] <- NULL
      if (!is.null(MAEQ.cv)) MAEQ.cv[[t]] <- NULL
      
      next
    }
    
    cID <- Last.time[, 1]
    
    val.cdata.pred <- val.cdata[
      match(cID, val.cdata[[ID]]),
      ,
      drop = FALSE
    ]
    
    keep <- !is.na(val.cdata.pred[[ID]])
    
    if (any(!keep)) {
      cID <- cID[keep]
      val.cdata.pred <- val.cdata.pred[keep, , drop = FALSE]
      survfit$Pred <- survfit$Pred[keep]
    }
    
    if (nrow(val.cdata.pred) == 0) {
      msg <- paste0(
        "The ", t,
        " th validation fold has no matched prediction subjects."
      )
      writeLines(msg)
      failed.folds[[as.character(t)]] <- msg
      
      if (!is.null(AUC.cv)) AUC.cv[[t]] <- NULL
      if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- NULL
      if (!is.null(Brier.cv)) Brier.cv[[t]] <- NULL
      if (!is.null(MAE.cv)) MAE.cv[[t]] <- NULL
      if (!is.null(MAEQ.cv)) MAEQ.cv[[t]] <- NULL
      
      next
    }
    
    if (CompetingRisk) {
      CIF <- as.data.frame(matrix(0, nrow = nrow(val.cdata.pred), ncol = 3))
      colnames(CIF) <- c("ID", "CIF1", "CIF2")
      
      CIF$ID <- val.cdata.pred[[ID]]
      CIF$time <- val.cdata.pred[, surv.var[1]]
      CIF$status <- val.cdata.pred[, surv.var[2]]
      
      if ("Brier Score" %in% metrics || "MAE" %in% metrics) {
        fitKM <- try(
          survival::survfit(New.surv.formula, data = val.cdata.pred),
          silent = TRUE
        )
        
        if ("try-error" %in% class(fitKM)) {
          msg <- paste0(
            "The censoring KM model failed in fold ",
            t,
            ". Brier/MAE will be set to NA for this fold."
          )
          writeLines(msg)
          fitKM <- NULL
          Gs <- NA
        } else {
          Gs.out <- try(
            summary(fitKM, times = landmark.time)$surv,
            silent = TRUE
          )
          
          if ("try-error" %in% class(Gs.out) || length(Gs.out) == 0) {
            Gs <- NA
          } else {
            Gs <- Gs.out
          }
        }
      }
      
      mean.AUC <- if ("AUC" %in% metrics) {
        matrix(NA, nrow = length(horizon.time), ncol = 2,
               dimnames = list(horizon.time, c("event1", "event2")))
      } else NULL
      
      mean.Cindex <- if ("Cindex" %in% metrics) {
        matrix(NA, nrow = length(horizon.time), ncol = 2,
               dimnames = list(horizon.time, c("event1", "event2")))
      } else NULL
      
      mean.Brier <- if ("Brier Score" %in% metrics) {
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
          ROC1 <- try(
            timeROC(
              T = CIF$time,
              delta = CIF$status,
              weighting = "marginal",
              marker = CIF$CIF1,
              cause = 1,
              times = horizon.time[j]
            ),
            silent = TRUE
          )
          
          if (!("try-error" %in% class(ROC1))) {
            mean.AUC[j, 1] <- ROC1$AUC_1[2]
          }
          
          CIF$status2 <- ifelse(
            CIF$status == 2,
            1,
            ifelse(CIF$status == 1, 2, 0)
          )
          
          ROC2 <- try(
            timeROC(
              T = CIF$time,
              delta = CIF$status2,
              weighting = "marginal",
              marker = CIF$CIF2,
              cause = 1,
              times = horizon.time[j]
            ),
            silent = TRUE
          )
          
          if (!("try-error" %in% class(ROC2))) {
            mean.AUC[j, 2] <- ROC2$AUC_1[2]
          }
        }
        
        if ("Cindex" %in% metrics) {
          cindex1 <- try(
            CindexCR(CIF$time, CIF$status, CIF$CIF1, 1),
            silent = TRUE
          )
          
          cindex2 <- try(
            CindexCR(CIF$time, CIF$status, CIF$CIF2, 2),
            silent = TRUE
          )
          
          if (!("try-error" %in% class(cindex1))) {
            mean.Cindex[j, 1] <- cindex1
          }
          
          if (!("try-error" %in% class(cindex2))) {
            mean.Cindex[j, 2] <- cindex2
          }
        }
        
        if ("Brier Score" %in% metrics || "MAE" %in% metrics) {
          if (is.null(fitKM) || is.na(Gs)) {
            
            if ("Brier Score" %in% metrics) {
              mean.Brier[j, 1] <- NA
              mean.Brier[j, 2] <- NA
            }
            
            if ("MAE" %in% metrics) {
              mean.MAE[j, 1] <- NA
              mean.MAE[j, 2] <- NA
            }
            
          } else {
            fitKM.horizon <- try(
              summary(fitKM, times = horizon.time[j]),
              silent = TRUE
            )
            
            if ("try-error" %in% class(fitKM.horizon) ||
                length(fitKM.horizon$surv) == 0) {
              
              if ("Brier Score" %in% metrics) {
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
                  val.cdata.pred[i, surv.var[1]] <= horizon.time[j] &&
                    val.cdata.pred[i, surv.var[2]] == 1
                )
                
                N2[i] <- as.numeric(
                  val.cdata.pred[i, surv.var[1]] <= horizon.time[j] &&
                    val.cdata.pred[i, surv.var[2]] == 2
                )
                
                if (val.cdata.pred[i, surv.var[1]] <= horizon.time[j] &&
                    val.cdata.pred[i, surv.var[2]] != 0) {
                  
                  Gt.out <- try(
                    summary(
                      fitKM,
                      times = as.numeric(val.cdata.pred[i, surv.var[1]])
                    )$surv,
                    silent = TRUE
                  )
                  
                  if ("try-error" %in% class(Gt.out) ||
                      length(Gt.out) == 0) {
                    Gt[i] <- NA
                  } else {
                    Gt[i] <- Gt.out
                  }
                  
                } else {
                  Gt[i] <- NA
                }
                
                if (val.cdata.pred[i, surv.var[1]] > horizon.time[j]) {
                  W.IPCW[i] <- 1 / (Gu / Gs)
                } else if (val.cdata.pred[i, surv.var[1]] <= horizon.time[j] &&
                           val.cdata.pred[i, surv.var[2]] != 0) {
                  W.IPCW[i] <- 1 / (Gt[i] / Gs)
                } else {
                  W.IPCW[i] <- NA
                }
              }
              
              RAWData.Brier <- data.frame(CIF, N1, N2, W.IPCW)
              
              if ("Brier Score" %in% metrics) {
                RAWData.Brier$Brier1 <- RAWData.Brier$W.IPCW *
                  abs(RAWData.Brier$CIF1 - RAWData.Brier$N1)^2
                
                RAWData.Brier$Brier2 <- RAWData.Brier$W.IPCW *
                  abs(RAWData.Brier$CIF2 - RAWData.Brier$N2)^2
                
                mean.Brier[j, 1] <- sum(
                  RAWData.Brier$Brier1,
                  na.rm = TRUE
                ) / nrow(RAWData.Brier)
                
                mean.Brier[j, 2] <- sum(
                  RAWData.Brier$Brier2,
                  na.rm = TRUE
                ) / nrow(RAWData.Brier)
              }
              
              if ("MAE" %in% metrics) {
                RAWData.Brier$MAE1 <- RAWData.Brier$W.IPCW *
                  abs(RAWData.Brier$CIF1 - RAWData.Brier$N1)
                
                RAWData.Brier$MAE2 <- RAWData.Brier$W.IPCW *
                  abs(RAWData.Brier$CIF2 - RAWData.Brier$N2)
                
                mean.MAE[j, 1] <- sum(
                  RAWData.Brier$MAE1,
                  na.rm = TRUE
                ) / nrow(RAWData.Brier)
                
                mean.MAE[j, 2] <- sum(
                  RAWData.Brier$MAE2,
                  na.rm = TRUE
                ) / nrow(RAWData.Brier)
              }
            }
          }
        }
        
        if ("MAEQ" %in% metrics) {
          
          quant1 <- quantile(
            CIF$CIF1,
            probs = seq(0, 1, by = quantile.width),
            na.rm = TRUE
          )
          
          EmpiricalCIF1 <- rep(NA, groups)
          PredictedCIF1 <- rep(NA, groups)
          
          for (i in seq_len(groups)) {
            
            subquant <- CIF[
              CIF$CIF1 > quant1[i] &
                CIF$CIF1 <= quant1[i + 1],
              1:2,
              drop = FALSE
            ]
            
            if (nrow(subquant) == 0) {
              EmpiricalCIF1[i] <- NA
              PredictedCIF1[i] <- NA
              next
            }
            
            quantsubdata <- val.cdata.pred[
              val.cdata.pred[[ID]] %in% subquant$ID,
              surv.var,
              drop = FALSE
            ]
            
            quantsubCIF <- try(
              GetEmpiricalCIF(
                data = quantsubdata,
                time = surv.var[1],
                status = surv.var[2]
              ),
              silent = TRUE
            )
            
            if ("try-error" %in% class(quantsubCIF)) {
              EmpiricalCIF1[i] <- NA
              PredictedCIF1[i] <- mean(subquant$CIF1, na.rm = TRUE)
              next
            }
            
            quantsubRisk1 <- quantsubCIF$H1
            ii <- 1
            
            while (ii <= nrow(quantsubRisk1)) {
              if (quantsubRisk1[ii, 1] > horizon.time[j]) {
                EmpiricalCIF1[i] <- if (ii >= 2) {
                  quantsubRisk1[ii - 1, 4]
                } else {
                  0
                }
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
            
            PredictedCIF1[i] <- mean(subquant$CIF1, na.rm = TRUE)
          }
          
          AllCIF1[[j]] <- data.frame(
            EmpiricalCIF1 = EmpiricalCIF1,
            PredictedCIF1 = PredictedCIF1
          )
          
          quant2 <- quantile(
            CIF$CIF2,
            probs = seq(0, 1, by = quantile.width),
            na.rm = TRUE
          )
          
          EmpiricalCIF2 <- rep(NA, groups)
          PredictedCIF2 <- rep(NA, groups)
          
          for (i in seq_len(groups)) {
            
            subquant <- CIF[
              CIF$CIF2 > quant2[i] &
                CIF$CIF2 <= quant2[i + 1],
              c(1, 3),
              drop = FALSE
            ]
            
            if (nrow(subquant) == 0) {
              EmpiricalCIF2[i] <- NA
              PredictedCIF2[i] <- NA
              next
            }
            
            quantsubdata <- val.cdata.pred[
              val.cdata.pred[[ID]] %in% subquant$ID,
              surv.var,
              drop = FALSE
            ]
            
            quantsubCIF <- try(
              GetEmpiricalCIF(
                data = quantsubdata,
                time = surv.var[1],
                status = surv.var[2]
              ),
              silent = TRUE
            )
            
            if ("try-error" %in% class(quantsubCIF)) {
              EmpiricalCIF2[i] <- NA
              PredictedCIF2[i] <- mean(subquant$CIF2, na.rm = TRUE)
              next
            }
            
            quantsubRisk2 <- quantsubCIF$H2
            ii <- 1
            
            while (ii <= nrow(quantsubRisk2)) {
              if (quantsubRisk2[ii, 1] > horizon.time[j]) {
                EmpiricalCIF2[i] <- if (ii >= 2) {
                  quantsubRisk2[ii - 1, 4]
                } else {
                  0
                }
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
            
            PredictedCIF2[i] <- mean(subquant$CIF2, na.rm = TRUE)
          }
          
          AllCIF2[[j]] <- data.frame(
            EmpiricalCIF2 = EmpiricalCIF2,
            PredictedCIF2 = PredictedCIF2
          )
        }
      }
      
      if (!is.null(AUC.cv)) AUC.cv[[t]] <- mean.AUC
      if (!is.null(Cindex.cv)) Cindex.cv[[t]] <- mean.Cindex
      if (!is.null(Brier.cv)) Brier.cv[[t]] <- mean.Brier
      if (!is.null(MAE.cv)) MAE.cv[[t]] <- mean.MAE
      
      if (!is.null(MAEQ.cv)) {
        names(AllCIF1) <- names(AllCIF2) <- horizon.time
        MAEQ.cv[[t]] <- list(
          AllCIF1 = AllCIF1,
          AllCIF2 = AllCIF2
        )
      }
    } else {
      Surv <- as.data.frame(matrix(0, nrow = nrow(val.cdata.pred), ncol = 2))
      colnames(Surv) <- c("ID", "Surv")
      Surv$ID <- val.cdata.pred[, ID]
      Surv$time <- val.cdata.pred[, surv.var[1]]
      Surv$status <- val.cdata.pred[, surv.var[2]]
      
      if ("Brier Score" %in% metrics || "MAE" %in% metrics) {
        fitKM <- try(
          survival::survfit(New.surv.formula, data = val.cdata.pred),
          silent = TRUE
        )
        
        if ("try-error" %in% class(fitKM)) {
          msg <- paste0(
            "The censoring KM model failed in fold ",
            t,
            ". Brier/MAE will be set to NA for this fold."
          )
          writeLines(msg)
          fitKM <- NULL
          Gs <- NA
        } else {
          Gs.out <- try(
            summary(fitKM, times = landmark.time)$surv,
            silent = TRUE
          )
          
          if ("try-error" %in% class(Gs.out) || length(Gs.out) == 0) {
            Gs <- NA
          } else {
            Gs <- Gs.out
          }
        }
      }
      
      mean.AUC <- if ("AUC" %in% metrics) {
        matrix(NA, nrow = length(horizon.time), ncol = 1,
               dimnames = list(horizon.time, "event"))
      } else NULL
      
      mean.Cindex <- if ("Cindex" %in% metrics) {
        matrix(NA, nrow = length(horizon.time), ncol = 1,
               dimnames = list(horizon.time, "event"))
      } else NULL
      
      mean.Brier <- if ("Brier Score" %in% metrics) {
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
        
        if ("Brier Score" %in% metrics || "MAE" %in% metrics) {
          
          fitKM.horizon <- try(summary(fitKM, times = horizon.time[j]), silent = TRUE)
          
          if ("try-error" %in% class(fitKM.horizon)) {
            
            if ("Brier Score" %in% metrics) mean.Brier[j, 1] <- NA
            if ("MAE" %in% metrics) mean.MAE[j, 1] <- NA
            
          } else {
            
            Gu <- fitKM.horizon$surv
            
            N1 <- vector()
            Gt <- vector()
            W.IPCW <- vector()
            
            for (i in seq_len(nrow(Surv))) {
              
              if (val.cdata.pred[i, surv.var[1]] <= horizon.time[j] &&
                  val.cdata.pred[i, surv.var[2]] == 1) {
                N1[i] <- 1
                Gt[i] <- summary(fitKM, times = val.cdata.pred[i, surv.var[1]])$surv
              } else {
                N1[i] <- 0
                Gt[i] <- NA
              }
              
              if (val.cdata.pred[i, surv.var[1]] > horizon.time[j]) {
                W.IPCW[i] <- 1 / (Gu / Gs)
              } else if (val.cdata.pred[i, surv.var[1]] <= horizon.time[j] &&
                         val.cdata.pred[i, surv.var[2]] == 1) {
                W.IPCW[i] <- 1 / (Gt[i] / Gs)
              } else {
                W.IPCW[i] <- NA
              }
            }
            
            RAWData.Brier <- data.frame(Surv, N1, W.IPCW)
            
            if ("Brier Score" %in% metrics) {
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
            
            quantsubdata <- val.cdata.pred[
              val.cdata.pred[, ID] %in% subquant$ID,
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
    jm.class = "mvjmcs",
    n.cv = n.cv,
    landmark.time = landmark.time,
    horizon.time = horizon.time,
    CompetingRisk = CompetingRisk,
    seed = seed,
    metrics = metrics,
    quantile.width = quantile.width,
    AUC.cv = AUC.cv,
    Cindex.cv = Cindex.cv,
    Brier.cv = Brier.cv,
    MAE.cv = MAE.cv,
    MAEQ.cv = MAEQ.cv,
    failed.folds = failed.folds
  )
  
  class(result) <- "DynPredAcc"
  
  return(result)
}