##' @title Time-dependent AUC  for joint models
##' @name AUCjmcs
##' @aliases AUCjmcs
##' @param seed a numeric value of seed to be specified for cross validation.
##' @param object object of class 'jmcs'.
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
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}, \link{survfitjmcs}}
##' @export
##' 

AUCjmcs <- function(seed = 100, object, landmark.time = NULL, horizon.time = NULL, 
                   obs.time = NULL, method = c("Laplace", "GH"), 
                   quadpoint = NULL, maxiter = NULL, n.cv = 3, 
                   survinitial = TRUE, ...) {
  
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
  CompetingRisk <- object$CompetingRisk
  set.seed(seed)
  cdata <- object$cdata
  ydata <- object$ydata
  long.formula <- object$LongitudinalSubmodel
  surv.formula <- object$SurvivalSubmodel
  surv.var <- all.vars(surv.formula)
  New.surv.formula.out <- paste0("survival::Surv(", surv.var[1], ",", 
                                 surv.var[2], "==0)")
  New.surv.formula <- as.formula(paste(New.surv.formula.out, 1, sep = "~"))
  random <- all.vars(object$random) 
  ID <- random[length(random)]
  
  folds <- caret::groupKFold(c(1:nrow(cdata)), k = n.cv)
  AUC.cv <- list()
  for (t in 1:n.cv) {
    
    train.cdata <- cdata[folds[[t]], ]
    train.ydata <- ydata[ydata[, ID] %in% train.cdata[, ID], ]
    
    fit <- try(jmcs(cdata = train.cdata, ydata = train.ydata, 
                    long.formula = long.formula,
                    surv.formula = surv.formula,
                    quadpoint = quadpoint, random = object$random, 
                    survinitial = survinitial), silent = TRUE)
    
    if ('try-error' %in% class(fit)) {
      writeLines(paste0("Error occured in the ", t, " th training!"))
      AUC.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      AUC.cv[[t]] <- NULL
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
        AUC.cv[[t]] <- NULL
      } else {
        if (CompetingRisk) {
          
          CIF <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 3))
          colnames(CIF) <- c("ID", "CIF1", "CIF2")
          CIF$ID <- val.cdata[, ID]
          CIF$time <- val.cdata[, surv.var[1]]
          CIF$status <- val.cdata[, surv.var[2]]
          mean.AUC <- matrix(NA, nrow = length(horizon.time), ncol = 2)
          for (j in 1:length(horizon.time)) {
            
              ## extract estimated CIF
              for (k in 1:nrow(CIF)) {
                CIF[k, 2] <- survfit$Pred[[k]][j, 2]
                CIF[k, 3] <- survfit$Pred[[k]][j, 3]
              }
              
            ROC <- timeROC::timeROC(T = CIF$time, delta = CIF$status,
                                    weighting = "marginal",
                                    marker = CIF$CIF1, cause = 1,
                                    times = horizon.time[j])
            
            mean.AUC[j, 1] <- ROC$AUC_1[2]
            
            CIF$status2 <- ifelse(CIF$status == 2, 1,
                                  ifelse(CIF$status == 1, 2, 0)
            )
            
            ROC <- timeROC::timeROC(T = CIF$time, delta = CIF$status2,
                                    weighting = "marginal",
                                    marker = CIF$CIF2, cause = 1,
                                    times = horizon.time[j])
            
            mean.AUC[j, 2] <- ROC$AUC_1[2]
            
          }
          
          AUC.cv[[t]] <- mean.AUC
            
          
        } else {
          
          Surv <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 2))
          colnames(Surv) <- c("ID", "Surv")
          Surv$ID <- val.cdata[, ID]
          Surv$time <- val.cdata[, surv.var[1]]
          Surv$status <- val.cdata[, surv.var[2]]
          mean.AUC <- matrix(NA, nrow = length(horizon.time), ncol = 1)
          for (j in 1:length(horizon.time)) {
            
            ## extract estimated Survival probability
            for (k in 1:nrow(Surv)) {
              Surv[k, 2] <- survfit$Pred[[k]][j, 2]
            }
            
            ROC <- timeROC::timeROC(T = Surv$time, delta = Surv$status,
                                    weighting = "marginal",
                                    marker = -Surv$Surv, cause = 1,
                                    times = horizon.time[j])
            
            mean.AUC[j, 1] <- ROC$AUC[2]
            
          }
          AUC.cv[[t]] <- mean.AUC
          
        }
        
      }
      
    }
    writeLines(paste0("The ", t, " th validation is done!"))
    
  }
  result <- list(n.cv = n.cv, AUC.cv = AUC.cv, landmark.time = landmark.time,
                 horizon.time = horizon.time, method = method, quadpoint = quadpoint, 
                 CompetingRisk = CompetingRisk, seed = seed)
  class(result) <- "AUCjmcs"
  
  return(result)
  
}