Concordance.mvjmcs <- function(seed = 100, object, n.cv = 3, maxiter = 10000,
                               initial.para = TRUE, ...) {
  
  if (!inherits(object, "mvjmcs"))
    stop("Use only with 'mvjmcs' xs.\n")
  
  CompetingRisk <- object$CompetingRisk
  set.seed(seed)
  cdata <- object$cdata
  ydata <- object$ydata
  long.formula <- object$LongitudinalSubmodel
  surv.formula <- object$SurvivalSubmodel
  surv.var <- all.vars(surv.formula)
  random <- object$random
  
  # ---- Longitudinal setup ----
  if(is.list(long.formula)){
    numBio = length(long.formula)
  } else {
    numBio = 1
  }
  random.form <- RE <- model <- vector("list", numBio)
  pREvec <- c()
  if(is.list(random)){
    for(g in 1:numBio){
      random.form[[g]] <- all.vars(random[[g]])
      pREvec[g] <- length(random.form[[g]])
      if(length(random.form[[g]])==1){
        # RE[[g]] <- NULL # already null
        model[[g]] <- "intercept"
      } else {
        RE[[g]] <- random.form[[g]][-length(all.vars(random[[g]]))]
        model[[g]] <- "interslope"
      }
    }
    ID <- random.form[[g]][length(all.vars(random[[g]]))]
  } else {
    for(g in 1:numBio){
      random.form[[g]] <- all.vars(random)
      pREvec[g] <- length(random.form[[g]])
      if (length(random.form[[g]]) == 1) {
        model[[g]] <- "intercept"
      } else {
        RE[[g]] <- random.form[[g]][-length(random.form[[g]])]
        model[[g]] <- "interslope"
      }
    }
    ID <- random.form[[g]][length(all.vars(random))]
  }
  
  if (initial.para & CompetingRisk) {
    beta <- object$beta
    initial.beta <- list()
    beta_names <- names(beta)
    bio_dim <- table(sub(".*_(bio\\d+)$", "\\1", beta_names))
    Start <- 0
    for (g in 1:numBio) {
      initial.beta[[g]] <- beta[(Start+1):(Start+bio_dim[g])]
      Start <- Start + bio_dim[g]
    }
    initial.alpha1 = object$alpha1
    initial.alpha2 = object$alpha2
    initial.alpha <- list(list(), list())
    Start <- 0
    for (g in 1:numBio) {
      initial.alpha[[1]][[g]] <- initial.alpha1[(Start+1):(Start+pREvec[g])]
      initial.alpha[[2]][[g]] <- initial.alpha2[(Start+1):(Start+pREvec[g])]
      Start <- Start + pREvec[g]
    }
    initial.para <- list(beta = initial.beta,
                         sigma = object$sigma, 
                         gamma1 = object$gamma1,
                         gamma2 = object$gamma2,
                         alpha = initial.alpha,
                         Sig = object$Sig)
  } else if (initial.para & !CompetingRisk) {
    beta <- object$beta
    initial.beta <- list()
    beta_names <- names(beta)
    bio_dim <- table(sub(".*_(bio\\d+)$", "\\1", beta_names))
    Start <- 0
    for (g in 1:numBio) {
      initial.beta[[g]] <- beta[(Start+1):(Start+bio_dim[g])]
      Start <- Start + bio_dim[g]
    }
    initial.alpha1 = object$alpha1
    initial.alpha <- list(list())
    Start <- 0
    for (g in 1:numBio) {
      initial.alpha[[1]][[g]] <- initial.alpha1[(Start+1):(Start+pREvec[g])]
      Start <- Start + pREvec[g]
    }
    initial.para <- list(beta = initial.beta,
                         sigma = object$sigma, 
                         gamma1 = object$gamma1,
                         alpha = initial.alpha,
                         Sig = object$Sig)
  } else {
    initial.para <- NULL
  }
  
  folds <- caret::groupKFold(c(1:nrow(cdata)), k = n.cv)
  Concordance.cv <- PI.cv <- list()
  for (t in 1:n.cv) {
    
    train.cdata <- cdata[folds[[t]], ]
    train.ydata <- ydata[ydata[, ID] %in% train.cdata[, ID], ]
    
    n.cores <- parallel::detectCores()
    fit <- try(
      mvjmcs(
        cdata = train.cdata,
        ydata = train.ydata,
        long.formula = long.formula,
        surv.formula = surv.formula,
        random = random,
        control = mvjmcs_control(
          maxiter = maxiter,
          tol = object$tol,
          initial.para = initial.para,
          opt = "optim",
          cpu.cores = n.cores
        )
      ),
      silent = TRUE
    )
    
    if ('try-error' %in% class(fit)) {
      writeLines(paste0("Error occured in the ", t, " th training!"))
      Concordance.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      Concordance.cv[[t]] <- NULL
    } else {
      
      val.cdata <- cdata[-folds[[t]], ]
      val.ydata <- ydata[ydata[, ID] %in% val.cdata[, ID], ]
      
      getinit <- Getmvinit(cdata = val.cdata, ydata = val.ydata, long.formula = long.formula,
                           surv.formula = surv.formula,
                           model = model, ID = ID, RE = RE,
                           REML = TRUE, random = random, opt = opt, initial.para, latAsso = "sre", landmark = FALSE)
      
      val.cdata2 <- getinit$cdata
      val.ydata2 <- getinit$ydata
      
      mdataM <- mdataSM <- vector("list", numBio)
      
      for(g in 1:numBio){
        subydata <- val.ydata2[[g]]
        submdata <- getinit$mdata[[g]]
        mdataM[[g]] <- submdata$ni
        n <- nrow(submdata)
        mdataSM[[g]] <- rep(0,n)
        submdata <- as.data.frame(submdata)
        submdata <- as.vector(submdata$ni)
        mdataSM[[g]][1] <- 1
        mdataCum <- cumsum(submdata)
        submdata2 <- submdata - 1
        mdataSM[[g]][2:n] <- mdataCum[2:n] - submdata2[2:n]
      }
      
      ## extract parameters
      
      if (CompetingRisk) {
        beta <- fit$betaList
        sigma <- fit$sigma
        gamma1 <- fit$gamma1
        gamma2 <- fit$gamma2
        alpha1 <- fit$alpha1
        alpha2 <- fit$alpha2
        alpha <- list(list(), list())
        Start <- 0
        for (g in 1:numBio) {
          alpha[[1]][[g]] <- alpha1[(Start+1):(Start+pREvec[g])]
          alpha[[2]][[g]] <- alpha2[(Start+1):(Start+pREvec[g])]
          Start <- Start + pREvec[g]
        }
        Sig <- fit$Sig
        H01 <- fit$H01
        H02 <- fit$H02
        
        Z <- getinit$Z
        X1 <- getinit$X1
        Y <- getinit$Y
        W <- getinit$W
        survtime <- getinit$survtime
        cmprsk <- getinit$cmprsk
        
        Postb <- GetBayes.mv(beta, sigma, gamma1, gamma2, alpha, H01, H02, 
                             Sig, Z, X1, Y, W, survtime, cmprsk, mdataM, mdataSM, "BFGS")
        
        X <- cbind(W, Postb)
        para1 <- c(gamma1, alpha1)
        para2 <- c(gamma2, alpha2)
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
        
        beta <- fit$betaList
        sigma <- fit$sigma
        gamma1 <- fit$gamma1
        alpha1 <- fit$alpha1
        alpha <- list(list())
        Start <- 0
        for (g in 1:numBio) {
          alpha[[1]][[g]] <- alpha1[(Start+1):(Start+pREvec[g])]
          Start <- Start + pREvec[g]
        }
        Sig <- fit$Sig
        H01 <- fit$H01
        
        Z <- getinit$Z
        X1 <- getinit$X1
        Y <- getinit$Y
        W <- getinit$W
        survtime <- getinit$survtime
        cmprsk <- getinit$cmprsk
        
        Postb <- GetBayesSF.mv(beta, sigma, gamma1, alpha, H01,
                                   Sig, Z, X1, Y, W, survtime, cmprsk, mdataM, mdataSM, "BFGS")
        
        X <- cbind(W, Postb)
        para1 <- c(gamma1, alpha1)
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
  class(result) <- "Concordance"
  
  return(result)
  
}