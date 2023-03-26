##' @title A metric of prediction accuracy of joint model by comparing the predicted risk
##' with the empirical risks stratified on different predicted risk group.
##' @name MAEQjmcs
##' @aliases MAEQjmcs
##' @param object object of class 'MAEQjmcs'.
##' @param seed a numeric value of seed to be specified for cross validation.
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
##' @param quantile.width a numeric value of width of quantile to be specified. Default is 0.25.
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}, \link{survfitjmcs}}
##' @export
##' 

MAEQjmcs <- function(object, seed = 100, landmark.time = NULL, horizon.time = NULL, 
                       obs.time = NULL, method = c("Laplace", "GH"), 
                       quadpoint = NULL, maxiter = 1000, 
                       n.cv = 3, survinitial = TRUE, 
                       quantile.width = 0.25, ...) {
  
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
  groups <- 1/quantile.width
  if (floor(groups) != groups)
    stop("The reciprocal of quantile.width must be an integer.")
  CompetingRisk <- object$CompetingRisk
  set.seed(seed)
  cdata <- object$cdata
  ydata <- object$ydata
  long.formula <- object$LongitudinalSubmodel
  surv.formula <- object$SurvivalSubmodel
  surv.var <- all.vars(surv.formula)
  random <- all.vars(object$random) 
  ID <- random[length(random)]
  
  folds <- caret::groupKFold(c(1:nrow(cdata)), k = n.cv)
  MAEQ.cv <- list()
  for (t in 1:n.cv) {
    
    train.cdata <- cdata[folds[[t]], ]
    train.ydata <- ydata[ydata[, ID] %in% train.cdata[, ID], ]
    
    fit <- try(jmcs(cdata = train.cdata, ydata = train.ydata, 
                      long.formula = long.formula,
                      surv.formula = surv.formula,
                      quadpoint = quadpoint, random = object$random,
                      maxiter = maxiter, 
                      survinitial = survinitial,
                    ), silent = TRUE)
    
    if ('try-error' %in% class(fit)) {
      writeLines(paste0("Error occured in the ", t, " th training!"))
      MAEQ.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      MAEQ.cv[[t]] <- NULL
    } else {
      
      val.cdata <- cdata[-folds[[t]], ]
      val.ydata <- ydata[ydata[, ID] %in% val.cdata[, ID], ]
      
      val.cdata <- val.cdata[val.cdata[, surv.var[1]] > landmark.time, ]
      val.ydata <- val.ydata[val.ydata[, ID] %in% val.cdata[, ID], ]
      val.ydata <- val.ydata[val.ydata[, obs.time] <= landmark.time, ]
      NewyID <- unique(val.ydata[, ID])
      val.cdata <- val.cdata[val.cdata[, ID] %in% NewyID, ]
      
      survfit <- try(survfitjmcs(fit, ynewdata = val.ydata, cnewdata = val.cdata, 
                                   u = horizon.time, method = method, 
                                   Last.time = rep(landmark.time, nrow(val.cdata)),
                                   obs.time = obs.time, quadpoint = quadpoint), silent = TRUE)
      
      if ('try-error' %in% class(survfit)) {
        writeLines(paste0("Error occured in the ", t, " th validation!"))
        MAEQ.cv[[t]] <- NULL
      } else {
        if (CompetingRisk) {
          AllCIF1 <- list()
          AllCIF2 <- list()
          for (j in 1:length(horizon.time)) {
            CIF <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 3))
            colnames(CIF) <- c("ID", "CIF1", "CIF2")
            CIF$ID <- val.cdata[, ID]
            ## extract estimated CIF
            for (k in 1:nrow(CIF)) {
              CIF[k, 2] <- survfit$Pred[[k]][j, 2]
              CIF[k, 3] <- survfit$Pred[[k]][j, 3]
            }
            ## group subjects based on CIF
            quant1 <- quantile(CIF$CIF1, probs = seq(0, 1, by = quantile.width))
            EmpiricalCIF1 <- rep(NA, groups)
            PredictedCIF1 <- rep(NA, groups)
            
            for (i in 1:groups) {
              subquant <- CIF[CIF$CIF1 > quant1[i] &
                                CIF$CIF1 <= quant1[i+1], 1:2]
              quantsubdata <- val.cdata[val.cdata[, ID] %in% subquant$ID, surv.var]
              
              quantsubCIF <- GetEmpiricalCIF(data = quantsubdata, 
                                             time = surv.var[1],
                                             status = surv.var[2])
              
              quantsubRisk1 <- quantsubCIF$H1
              ii <- 1
              while (ii <= nrow(quantsubRisk1)) {
                if (quantsubRisk1[ii, 1] > horizon.time[j]) {
                  if (ii >= 2) {
                    EmpiricalCIF1[i] <- quantsubRisk1[ii-1, 4]
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
                                CIF$CIF2 <= quant2[i+1], c(1, 3)]
              quantsubdata <- cdata[cdata[, ID] %in% subquant$ID, surv.var]
              
              quantsubCIF <- GetEmpiricalCIF(data = quantsubdata, 
                                             time = surv.var[1],
                                             status = surv.var[2])
              
              quantsubRisk2 <- quantsubCIF$H2
              ii <- 1
              while (ii <= nrow(quantsubRisk2)) {
                if (quantsubRisk2[ii, 1] > horizon.time[j]) {
                  if (ii >= 2) {
                    EmpiricalCIF2[i] <- quantsubRisk2[ii-1, 4]
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
          names(AllCIF1) <- names(AllCIF2) <- horizon.time
          result <- list(AllCIF1 = AllCIF1, AllCIF2 = AllCIF2)
        } else {
          AllSurv <- list()
          for (j in 1:length(horizon.time)) {
            Surv <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 2))
            colnames(Surv) <- c("ID", "Surv")
            Surv$ID <- val.cdata[, ID]
            ## extract estimated survival prob
            for (k in 1:nrow(Surv)) {
              Surv[k, 2] <- survfit$Pred[[k]][j, 2]
            }
            ## group subjects based on survival prob
            quant <- quantile(Surv$Surv, probs = seq(0, 1, by = quantile.width))
            EmpiricalSurv <- rep(NA, groups)
            PredictedSurv <- rep(NA, groups)
            for (i in 1:groups) {
              subquant <- Surv[Surv$Surv > quant[i] &
                                 Surv$Surv <= quant[i+1], c(1, 2)]
              quantsubdata <- cdata[cdata[, ID] %in% subquant$ID, surv.var]
              colnames(quantsubdata) <- c("time", "status")
              fitKM <- survfit(Surv(time, status) ~ 1, data = quantsubdata)
              fitKM.horizon <- try(summary(fitKM, times = horizon.time[j]), silent = TRUE)
              if ('try-error' %in% class(fitKM.horizon)) {
                EmpiricalSurv[i] <- summary(fitKM, times = max(quantsubdata$time))$surv
              } else {
                EmpiricalSurv[i] <- summary(fitKM, times = horizon.time[j])$surv
              }
              PredictedSurv[i] <-mean(subquant$Surv)
            }
            AllSurv[[j]] <- data.frame(EmpiricalSurv, PredictedSurv)
          }
          names(AllSurv) <- horizon.time
          result <- list(AllSurv = AllSurv)
        }
        MAEQ.cv[[t]] <- result
        writeLines(paste0("The ", t, " th validation is done!"))
      }
    }
  }
  result <- list(MAEQ.cv = MAEQ.cv, n.cv = n.cv, landmark.time = landmark.time,
                 horizon.time = horizon.time, method = method, quadpoint = quadpoint, 
                 CompetingRisk = CompetingRisk, seed = seed)
  class(result) <- "MAEQjmcs"
  return(result)
}
