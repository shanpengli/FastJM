##' @title Time-dependent AUC  for joint models
##' @name PAjmcs
##' @aliases PAjmcs
##' @param seed a numeric value of seed to be specified for cross validation.
##' @param object object of class 'jmcs'.
##' @param landmark.time a numeric value of time for which dynamic prediction starts..
##' @param horizon.time a numeric vector of future times for which predicted probabilities are to be computed.
##' @param obs.time a character string of specifying a longitudinal time variable.
##' @param quadpoint the number of pseudo-adaptive Gauss-Hermite quadrature points to be used.
##' @param maxiter the maximum number of iterations of the EM algorithm that the 
##' function will perform. Default is 10000.
##' @param n.cv number of folds for cross validation. Default is 3.
##' @param survinitial Fit a Cox model to obtain initial values of the parameter estimates. Default is TRUE.
##' @param initial.para Initial guess of parameters for cross validation. Default is FALSE.
##' @param LOCF a logical value to indicate whether the last-observation-carried-forward approach applies to prediction. 
##' If \code{TRUE}, then \code{LOCFcovariate} and \code{clongdata} must be specified to indicate 
##' which time-dependent survival covariates are included for dynamic prediction. Default is FALSE.
##' @param LOCFcovariate a vector of string with time-dependent survival covariates if \code{LOCF = TRUE}. Default is NULL.
##' @param clongdata a long format data frame where time-dependent survival covariates are incorporated. Default is NULL.
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}, \link{survfitjmcs}}
##' @export
##' 

PAjmcs <- function(seed = 100, object, landmark.time = NULL, horizon.time = NULL, 
                    obs.time = NULL, quadpoint = NULL, maxiter = NULL, n.cv = 3, 
                    survinitial = TRUE, initial.para = FALSE,
                    LOCF = FALSE, LOCFcovariate = NULL, clongdata = NULL, ...) {
  
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' xs.\n")
  if (is.null(landmark.time)) 
    stop("Please specify the landmark.time for dynamic prediction.")   
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
  
  folds <- caret::groupKFold(c(1:nrow(cdata)), k = n.cv)
  PA.cv <- list()
  for (t in 1:n.cv) {
    
    train.cdata <- cdata[folds[[t]], ]
    train.ydata <- ydata[ydata[, ID] %in% train.cdata[, ID], ]
    
    fit <- try(jmcs(cdata = train.cdata, ydata = train.ydata, 
                    long.formula = long.formula,
                    surv.formula = surv.formula,
                    quadpoint = quadpoint, random = object$random, 
                    survinitial = survinitial,
                    opt = object$opt, maxiter = maxiter,
                    initial.para = initial.para), silent = TRUE)
    
    if ('try-error' %in% class(fit)) {
      writeLines(paste0("Error occured in the ", t, " th training!"))
      PA.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      PA.cv[[t]] <- NULL
    } else {
      
      val.cdata <- cdata[-folds[[t]], ]
      val.ydata <- ydata[ydata[, ID] %in% val.cdata[, ID], ]
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
      
      survfit <- try(survPAjmcs(fit, ynewdata = val.ydata, cnewdata = val.cdata,
                                u = horizon.time, s = landmark.time,
                                obs.time = obs.time, quadpoint = quadpoint,
                                LOCF = LOCF, LOCFcovariate = LOCFcovariate, 
                                clongdata = val.clongdata), silent = TRUE)
      
      if ('try-error' %in% class(survfit)) {
        writeLines(paste0("Error occured in the ", t, " th validation!"))
        PA.cv[[t]] <- NULL
      } else {
        if (CompetingRisk) {
          
          mean.PA <- matrix(NA, nrow = length(horizon.time), ncol = 2)
          for (j in 1:length(horizon.time)) {
            mean.PA[j, 1] <- as.numeric(survfit$Pred$Risk1[[j]]$Psuedo.R)
            mean.PA[j, 2] <- as.numeric(survfit$Pred$Risk2[[j]]$Psuedo.R)
          }
          
          PA.cv[[t]] <- mean.PA
          
        } else {
          
          stop("TBD.")
          
        }
        
      }
      
    }
    writeLines(paste0("The ", t, " th validation is done!"))
    
  }
  result <- list(n.cv = n.cv, PA.cv = PA.cv, landmark.time = landmark.time,
                 horizon.time = horizon.time, quadpoint = quadpoint, 
                 CompetingRisk = CompetingRisk, seed = seed)
  class(result) <- "PAjmcs"
  
  return(result)
  
}