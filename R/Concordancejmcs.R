##' @title Concordance for joint models
##' @name Concordancejmcs
##' @aliases Concordancejmcs
##' @param seed a numeric value of seed to be specified for cross validation.
##' @param object object of class 'jmcs'.
##' @param opt Optimization method to fit a linear mixed effects model, either nlminb (default) or optim.
##' @param n.cv number of folds for cross validation. Default is 3.
##' @param maxiter the maximum number of iterations of the EM algorithm that the 
##' function will perform. Default is 10000.
##' @param initial.optimizer Method for numerical optimization to be used. Default is \code{BFGS}.
##' @param initial.para Initial guess of parameters for cross validation. Default is FALSE.
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}}
##' @export
##' 

Concordancejmcs <- function(seed = 100, object, opt = "nlminb", n.cv = 3, maxiter = 10000,
                              initial.optimizer = "BFGS", initial.para = TRUE, ...) {
  
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' xs.\n")
  
  CompetingRisk <- object$CompetingRisk
  set.seed(seed)
  cdata <- object$cdata
  ydata <- object$ydata
  long.formula <- object$LongitudinalSubmodel
  surv.formula <- object$SurvivalSubmodel
  surv.var <- all.vars(surv.formula)
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
  Concordance.cv <- PI.cv <- list()
  for (t in 1:n.cv) {
    
    train.cdata <- cdata[folds[[t]], ]
    train.ydata <- ydata[ydata[, ID] %in% train.cdata[, ID], ]
    
    fit <- try(jmcs(cdata = train.cdata, ydata = train.ydata, 
                      long.formula = long.formula,
                      surv.formula = surv.formula,
                      quadpoint = object$quadpoint, random = object$random, 
                      maxiter = maxiter, opt = opt, 
                      tol = object$tol, initial.para = initial.para), silent = TRUE)
    
    if ('try-error' %in% class(fit)) {
      writeLines(paste0("Error occured in the ", t, " th training!"))
      Concordance.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      Concordance.cv[[t]] <- NULL
    } else {
      
      val.cdata <- cdata[-folds[[t]], ]
      val.ydata <- ydata[ydata[, ID] %in% val.cdata[, ID], ]
      
      getinit <- Getinit(cdata = val.cdata, ydata = val.ydata, long.formula = long.formula,
                         surv.formula = surv.formula,
                         model = model, ID = ID, RE = RE, random = object$random, survinitial = FALSE,
                         initial.para = initial.para, opt = opt)
      
      val.cdata2 <- getinit$cdata
      val.ydata2 <- getinit$ydata
      
      ## extract parameters
      
      if (CompetingRisk) {
        beta <- fit$beta
        sigma <- fit$sigma
        gamma1 <- fit$gamma1
        gamma2 <- fit$gamma2
        nu1 <- fit$nu1
        nu2 <- fit$nu2
        Sig <- fit$Sig
        H01 <- fit$H01
        H02 <- fit$H02
        
        Z <- getinit$Z
        X1 <- getinit$X1
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
        
        Posttheta <- GetBayes(beta, sigma, gamma1, gamma2, nu1, nu2, H01, H02, 
                              Sig, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS, initial.optimizer)
        
        X <- cbind(X2, Posttheta)
        para1 <- c(gamma1, nu1)
        para2 <- c(gamma2, nu2)
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
      }
      
    }
    
  }
  result <- list(n.cv = n.cv, Concordance.cv = Concordance.cv, PI.cv = PI.cv,
                 CompetingRisk = CompetingRisk, seed = seed)
  class(result) <- "Concordancejmcs"
  
  return(result)
  
}