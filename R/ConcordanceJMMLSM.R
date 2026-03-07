##' @title Concordance for joint models
##' @name ConcordanceJMMLSM
##' @aliases ConcordanceJMMLSM
##' @param seed a numeric value of seed to be specified for cross validation.
##' @param object object of class 'JMMLSM'.
##' @param opt Optimization method to fit a linear mixed effects model, either nlminb (default) or optim.
##' @param n.cv number of folds for cross validation. Default is 3.
##' @param maxiter the maximum number of iterations of the EM algorithm that the 
##' function will perform. Default is 10000.
##' @param initial.optimizer Method for numerical optimization to be used. Default is \code{BFGS}.
##' @param initial.para Initial guess of parameters for cross validation. Default is FALSE.
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{JMMLSM}}
##' @export
##' 

ConcordanceJMMLSM <- function(seed = 100, object, opt = "nlminb", n.cv = 3, maxiter = 10000,
                              initial.optimizer = "BFGS", initial.para = TRUE, ...) {
  
  CompetingRisk <- object$CompetingRisk
  set.seed(seed)
  cdata <- object$cdata
  ydata <- object$ydata
  long.formula <- object$LongitudinalSubmodelmean
  surv.formula <- object$SurvivalSubmodel
  surv.var <- all.vars(surv.formula)
  variance.formula <- as.formula(paste("", object$LongitudinalSubmodelvariance[3], sep = "~"))
  random.form <- all.vars(object$random) 
  ID <- random.form[length(random.form)]
  
  if (length(random.form) == 1) {
    RE <- NULL
    model <- "intercept"
  } else {
    RE <- random.form[-length(random.form)]
    model <- "interslope"
  }
  
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
  Concordance.cv <- PI.cv <- list()
  for (t in 1:n.cv) {
    
    train.cdata <- cdata[folds[[t]], ]
    train.ydata <- ydata[ydata[, ID] %in% train.cdata[, ID], ]
    
    fit <- try(JMMLSM(cdata = train.cdata, ydata = train.ydata, 
                      long.formula = long.formula,
                      surv.formula = surv.formula,
                      variance.formula = variance.formula, 
                      quadpoint = object$quadpoint, random = object$random, 
                      maxiter = maxiter, opt = opt, 
                      epsilon = object$epsilon, initial.para = initial.para), silent = TRUE)
    
    if ('try-error' %in% class(fit)) {
      writeLines(paste0("Error occured in the ", t, " th training!"))
      Concordance.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      Concordance.cv[[t]] <- NULL
    } else {
      
      val.cdata <- cdata[-folds[[t]], ]
      val.ydata <- ydata[ydata[, ID] %in% val.cdata[, ID], ]
      
      getinit <- Getinit.JMH(cdata = val.cdata, ydata = val.ydata, long.formula = long.formula,
                         surv.formula = surv.formula, variance.formula = variance.formula,
                         model = model, ID = ID, RE = RE, random = object$random, survinitial = FALSE,
                         initial.para = initial.para, opt = opt)
      
      val.cdata2 <- getinit$cdata
      val.ydata2 <- getinit$ydata
    
      ## extract parameters
      
      if (CompetingRisk) {
        beta <- fit$beta
        tau <- fit$tau
        gamma1 <- fit$gamma1
        gamma2 <- fit$gamma2
        alpha1 <- fit$alpha1
        alpha2 <- fit$alpha2
        vee1 <- fit$vee1
        vee2 <- fit$vee2
        Sig <- fit$Sig
        H01 <- fit$H01
        H02 <- fit$H02
        
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
        
        Posttheta <- GetBayes.JMH(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, 
                             Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, initial.optimizer)
        
        X <- cbind(X2, Posttheta)
        para1 <- c(gamma1, alpha1, vee1)
        para2 <- c(gamma2, alpha2, vee2)
        Risk1Score <- X %*% para1
        Risk2Score <- X %*% para2
        subcdata <- data.frame(survtime, cmprsk, exp(Risk1Score), exp(Risk2Score))
        colnames(subcdata) <- c("time", "status", "Risk1Score", "Risk2Score")
        PI.cv[[t]] <- subcdata
        
        Risk1Cindex <- CindexCR(subcdata$time, subcdata$status, 
                                   subcdata$Risk1Score, Cause_int = 1)
        
        Risk2Cindex <- CindexCR(subcdata$time, subcdata$status, 
                                   subcdata$Risk2Score, Cause_int = 2)
        
        Concordance.cv[[t]] <- c(Risk1Cindex, Risk2Cindex)
        writeLines(paste0("The ", t, " th validation is done!"))
      } else {
        
        beta <- fit$beta
        tau <- fit$tau
        gamma1 <- fit$gamma1
        alpha1 <- fit$alpha1
        vee1 <- fit$vee1
        Sig <- fit$Sig
        H01 <- fit$H01
        
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
        
        Posttheta <- GetBayesSF.JMH(beta, tau, gamma1, alpha1, vee1, H01,
                                  Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, initial.optimizer)
        
        X <- cbind(X2, Posttheta)
        para1 <- c(gamma1, alpha1, vee1)
        Risk1Score <- X %*% para1
        subcdata <- data.frame(survtime, cmprsk, exp(Risk1Score))
        colnames(subcdata) <- c("time", "status", "Risk1Score")
        PI.cv[[t]] <- subcdata
        
        Risk1Cindex <- CindexCR(subcdata$time, subcdata$status, 
                                subcdata$Risk1Score, Cause_int = 1)
        
        Concordance.cv[[t]] <- Risk1Cindex
        writeLines(paste0("The ", t, " th validation is done!"))
        
      }
      
    }
    
  }
  result <- list(n.cv = n.cv, Concordance.cv = Concordance.cv, PI.cv = PI.cv,
                 CompetingRisk = CompetingRisk, seed = seed)
  class(result) <- "ConcordanceJMMLSM"
  
  return(result)
  
}