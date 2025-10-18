##' Joint modeling of multivariate longitudinal continuous data and competing risks
##'
##' Function fits a joint model for multiple longitudinal outcomes and competing risks using a fast EM algorithm.
##'
##' @title Joint modeling of multivariate longitudinal and competing risks data
##' @name mvjmcs
##' @param ydata A longitudinal data frame in long format.
##' @param cdata A survival data frame with competing risks or single failure. Each subject has one data entry.
##' @param long.formula A list of formula objects specifying fixed effects for each longitudinal outcome.
##' @param random A formula or list of formulas describing random effects structures (e.g., \code{~ 1|ID}).
##' @param surv.formula A formula for the survival sub-model, including survival time and event indicator.
##' @param maxiter Maximum number of EM iterations. Default is 10000.
##' @param opt Optimization method for mixed model. Default is \code{"nlminb"}.
##' @param tol Convergence tolerance for EM algorithm. Default is 0.0001.
##' @param print.para Logical; if \code{TRUE}, prints parameter values at each iteration.
##' @param initial.para Optional list of initialized parameters. Default is \code{NULL}.
##' @param cpu.cores Number of CPU cores for parallel computation. Default is 1. 
##'
##' @return  Object of class \code{mvjmcs} with elements
##' \item{beta}{the vector of all biomarker-specific fixed effects for the linear mixed effects sub-models.} 
##' \item{betaList}{the list of biomarker-specific fixed effects for the linear mixed effects sub-model.} 
##' \item{gamma1}{the vector of fixed effects for type 1 failure for the survival model.}
##' \item{gamma2}{the vector of fixed effects for type 2 failure for the survival model. 
##' Valid only if \code{CompetingRisk = TRUE}.}
##' \item{alpha1}{the vector of association parameter(s) for type 1 failure.}
##' \item{alpha2}{the vector of association parameter(s) for type 2 failure. Valid only if \code{CompetingRisk = TRUE}.}
##' \item{H01}{the matrix that collects baseline hazards evaluated at each uncensored event time for type 1 failure. 
##' The first column denotes uncensored event times, the second column the number of events, and the third columns 
##' the hazards obtained by Breslow estimator.}
##' \item{H02}{the matrix that collects baseline hazards evaluated at each uncensored event time for type 2 failure. 
##' The data structure is the same as \code{H01}. Valid only if \code{CompetingRisk = TRUE}.}
##' \item{Sig}{the variance-covariance matrix of the random effects.}
##' \item{sigma}{the vector of the variance of the biomarker-specific measurement error for the linear mixed effects sub-models.}
##' \item{iter}{the total number of iterations until convergence.}
##' \item{convergence}{convergence identifier: 1 corresponds to successful convergence, 
##' whereas 0 to a problem (i.e., when 0, usually more iterations are required).}
##' \item{vcov}{the variance-covariance matrix of all the fixed effects for both models.}
##' \item{FisherInfo}{the Empirical Fisher information matrix.}
##' \item{Score}{a matrix of the score function for all subjects.}
##' \item{sebeta}{the standard error of \code{beta}.}
##' \item{segamma1}{the standard error of \code{gamma1}.}
##' \item{segamma2}{the standard error of \code{gamma2}. 
##' Valid only if \code{CompetingRisk = TRUE}.}
##' \item{sealpha1}{the standard error of \code{nu1}.}
##' \item{sealpha2}{the standard error of \code{nu2}. Valid only if \code{CompetingRisk = TRUE}.}
##' \item{seSig}{the vector of standard errors of covariance of random effects.}
##' \item{sesigma}{the standard error of variance of biomarker-specific measurement error for the linear mixed effects sub-models.}
##' \item{pos.mode}{the posterior mode of the conditional distribution of random effects.}
##' \item{pos.cov}{the posterior covariance of the conditional distribution of random effects.}
##' \item{CompetingRisk}{logical value; TRUE if a competing event are accounted for.}
##' \item{ydata}{the input longitudinal dataset for fitting a joint model.
##' It has been re-ordered in accordance with descending observation times in \code{cdata}.}
##' \item{cdata}{the input survival dataset for fitting a joint model.
##' It has been re-ordered in accordance with descending observation times.}
##' \item{PropEventType}{a frequency table of number of events.}
##' \item{LongitudinalSubmodel}{the component of the \code{long.formula}.}
##' \item{SurvivalSubmodel}{the component of the \code{surv.formula}.}
##' \item{random}{the component of the \code{random}.}
##' \item{call}{the matched call.}
##' \item{id}{the grouping vector for the longitudinal outcome.}
##' \item{opt}{the numerical optimizer for obtaining the initial guess of the parameters in the linear mixed effects sub-models.}
##' \item{runtime}{the total computation time.}
##'
##' @examples
##' 
##' 
##'   require(FastJM)
##'   require(survival)
##'   require(future)
##'   require(future.apply)
##'   
##'   data(mvcdata)
##'   data(mvydata)
##' 
##'   \donttest{
##'   # Fit joint model with two biomarkers
##'   fit <- mvjmcs(ydata = mvydata, cdata = mvcdata, 
##'                 long.formula = list(Y1 ~ X11 + X12 + time, 
##'                                     Y2 ~ X11 + X12 + time),
##'                 random = list(~ time | ID,
##'                               ~ 1 | ID),
##'                 surv.formula = Surv(survtime, cmprsk) ~ X21 + X22, maxiter = 1000, opt = "optim", 
##'                 tol = 1e-3, print.para = FALSE)
##'   fit
##'   
##'   # Extract the parameter estimates of longitudinal sub-model fixed effects
##'   fixef(fit, process = "Longitudinal")
##'   
##'   # Extract the parameter estimates of survival sub-model fixed effects
##'   fixef(fit, process = "Event")
##'   
##'   # Obtain the random effects estimates for first 6 subjects 
##'   head(ranef(fit))
##'   }
##'   
##' @export

mvjmcs <- function(ydata, cdata, long.formula,
                   random = NULL, surv.formula,
                   maxiter = 10000, opt = "nlminb", tol = 0.005,
                   print.para = TRUE, 
                   initial.para = NULL, cpu.cores = NULL){
  
  start_time <- Sys.time()
  
  # ---- Longitudinal setup ----
  if(is.list(long.formula)){
    numBio = length(long.formula)
  } else {
    numBio = 1
  }
  
  random.form <- RE  <- model <- vector("list", numBio)
  
  if(is.list(random)){
    for(g in 1:numBio){
      random.form[[g]] <- all.vars(random[[g]])
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
      if(!is.list(random)){
        # RE[[g]] <- NULL
        model[[g]] <- "intercept"
      } else {
        RE[[g]] <- random.form[[g]][-length(all.vars(random))]
        model[[g]] <- "interslope"
      }
    }
    ID <- random.form[[g]][length(all.vars(random))]
  }
  
  lengthb <- length(long.formula)
  long <- list()
  
  # ---- Formula checks for each biomarker ----
  for(g in 1:lengthb){
    if (!inherits(long.formula[[g]], "formula") || length(long.formula[[g]]) != 3) {
      stop("\nLinear mixed effects model must be a formula of the form \"resp ~ pred\"")
    }
    long[[g]] <- all.vars(long.formula[[g]])
  }
  longfmla <- list(lengthb)
  for(g in 1:lengthb){
    longfmla[g] <- long.formula[g]
  }
  
  
  # ---- CPU setup ----
  if(is.null(cpu.cores)){
    cpu.cores <- 1
  }

  cnames <- colnames(cdata)
  ynames <- colnames(ydata)
  
  survfmla <- surv.formula
  
  rawydata <- ydata
  rawcdata <- cdata
  
  getinit <- Getmvinit(cdata = cdata, ydata = ydata, long.formula = long.formula,
                       surv.formula = surv.formula,
                       model = model, ID = ID, RE = RE,
                       REML = TRUE, random = random, opt = opt, initial.para)
  
  if (is.null(getinit)) {
    stop("Numerical failure occurred when fitting a linear mixed effects model for initial guess.")
  }
  
  
  mdataM <- mdataSM <- vector("list", numBio)
  
  for(g in 1:numBio){
    ydata <- getinit$ydata[[g]]
    mdata <- getinit$mdata[[g]]
    mdataM[[g]] <- mdata$ni
    n <- nrow(mdata)
    mdataSM[[g]] <- rep(0,n)
    mdata <- as.data.frame(mdata)
    mdata <- as.vector(mdata$ni)
    mdataSM[[g]][1] <- 1
    mdataCum <- cumsum(mdata)
    mdata2 <- mdata - 1
    mdataSM[[g]][2:n] <- mdataCum[2:n] - mdata2[2:n]
  }
  
  # mdataM <- getinit$mdata
  
  
  cdata <- getinit$cdata
  
  survival <- all.vars(surv.formula)
  status <- as.vector(cdata[, survival[2]])
  if (prod(c(0, 1, 2) %in% unique(status))) {
    ## initialize parameters
    
    getriskset <- Getriskset(cdata = cdata, surv.formula = surv.formula)
    
    ## number of distinct survival time
    H01 <- getriskset$tablerisk1
    H02 <- getriskset$tablerisk2
    
    CompetingRisk <- TRUE
  } else {
    getriskset <- Getriskset(cdata = cdata, surv.formula = surv.formula)
    
    ## number of distinct survival time
    H01 <- getriskset$tablerisk1
    CompetingRisk <- FALSE
  }
  
  # GH.val  <- gauss.quad.prob(quadpoint)
  # weight.c <- GH.val$weights # weights
  # abscissas.c <- GH.val$nodes # abscissas
  
  
  survtime <- getinit$survtime
  cmprsk <- getinit$cmprsk
  
  iter=0
  
  # ======================== #
  # ~~~~ Competing Risk ~~~~ #
  # ======================== #
  if(CompetingRisk == TRUE){
    
    pREvec <- c()
    
    for(g in 1:numBio){
      pREvec[g] <- ncol(getinit$Z[[g]])
    }
    pREtotal <- sum(pREvec)
    
    index = 0
    if (is.null(initial.para)) {
      
      SigList <- getinit$Sig
      Sig <- matrix(rep(0, pREtotal), nrow = pREtotal, ncol = pREtotal)
      for(g in 1:length(pREvec)){
        Sig[(index+1):(index+pREvec[g]),(index+1):(index+pREvec[g])] <- SigList[[g]]
        index = index + pREvec[g]
      }
      
    } else {
      Sig <- getinit$Sig
    }
    
    
    numSubj <- n
    # opt <- list()
    pos.mode <- vector("list", numSubj)
    subX1 <- subY <- subZ <- vector("list", numSubj)
    pos.cov <- list()
    namesbeta <- vector("list", numBio)
    pbeta <- c()
    for (g in 1:numBio) {
      namesbeta[[g]] <- paste0(names(getinit$beta[[g]]), "_bio", g)
      pbeta[g] <- length(getinit$beta[[g]])
    }
    namesgamma <- names(getinit$gamma1)
    
    beta <- unlist(getinit$beta)
    betaList <- getinit$beta
    gamma1 <- getinit$gamma1
    gamma2 <- getinit$gamma2
    alphaList <- getinit$alpha
    alpha1 <- unlist(alphaList[[1]])
    alpha2 <- unlist(alphaList[[2]])
    sigma <- getinit$sigma
    
    for(j in 1:numSubj){
      subX1[[j]] <- vector("list", numBio)
      for (g in 1:numBio) {
        
        numRep <- mdataM[[g]][j]
        indexStart <- mdataSM[[g]][j]
        
        subX1[[j]][[g]] <- getinit$X1[[g]][indexStart:(indexStart+numRep-1),, drop = FALSE]
        subY[[j]][[g]] <- getinit$Y[[g]][indexStart:(indexStart+numRep-1)]
        subZ[[j]][[g]] <- getinit$Z[[g]][indexStart:(indexStart+numRep-1),, drop = FALSE]
      }
    }
    
    survtime <- getinit$survtime
    cmprsk <- getinit$cmprsk
    W <- getinit$W
    ngamma <- ncol(W)
    
    old_parallelly <- options(
      parallelly.makeNodePSOCK.rscript_args = c("--vanilla")
    )
    on.exit(options(old_parallelly), add = TRUE)
    old_plan <- future::plan()
    future::plan(future::multisession, workers = cpu.cores)
    on.exit(future::plan(old_plan), add = TRUE)
    
    repeat{
      iter <- iter + 1
      
      CUH01 <- rep(0, n) 
      CUH02 <- rep(0, n) 
      HAZ01 <- rep(0, n) 
      HAZ02 <- rep(0, n) 
      CumuH01 <- cumsum(H01[, 3]) 
      CumuH02 <- cumsum(H02[, 3]) 
      getHazard(CumuH01, CumuH02, survtime, cmprsk, H01, H02, CUH01, CUH02, HAZ01, HAZ02) 
      
      prebeta <- beta
      pregamma1 <- gamma1
      pregamma2 <- gamma2
      prealphaList <- alphaList
      prealpha1 <- alpha1
      prealpha2 <- alpha2
      preH01 <- H01
      preH02 <- H02
      preSig <- Sig
      presigma <- sigma
      
      if (print.para) {
        writeLines("iter is:")
        print(iter)
        writeLines("beta is:")
        print(prebeta)
        writeLines("sigma is:")
        print(presigma)
        writeLines("Sig is:")
        print(preSig)
        writeLines("gamma1 is:")
        print(pregamma1)
        writeLines("gamma2 is:")
        print(pregamma2)
        writeLines("alpha1 is:")
        print(prealpha1)
        writeLines("alpha2 is:")
        print(prealpha2)
      }
      
      data <- list(
        beta = betaList, gamma1 = gamma1, gamma2 = gamma2,
        alpha = alphaList, sigma = sigma,
        Z = subZ, X1 = subX1, Y = subY, Sig = Sig,
        CUH01 = CUH01, CUH02 = CUH02, HAZ01 = HAZ01, HAZ02 = HAZ02,
        cmprsk = cmprsk, W = W
      )
      
      res <- future.apply::future_lapply(seq_len(numSubj), estepMV_worker, future.seed = TRUE, future.scheduling = 2,
                                         data, pREtotal)
      pos.mode <- lapply(res, `[[`, "mode")
      pos.cov  <- lapply(res, function(x) crossprod(x$ccov))
      
      output <- normalApprox(
        subX1, subY, subZ, W,
        mdataM, mdataSM,
        pos.mode, sigma, pos.cov,
        H01, H02, survtime, cmprsk,
        gamma1, gamma2, alphaList,
        CUH01, CUH02, HAZ01, HAZ02, Sig, betaList
      )
      
      beta <- output$beta
      betaList <- output$betaList
      sigma <- output$sigmaVec
      Sig <- output$Sig
      
      gamma1 <- output$phi1[1:ngamma]
      gamma2 <- output$phi2[1:ngamma]
      alpha1 <- output$phi1[(ngamma+1):(length(output$phi1))]
      alpha2 <- output$phi2[(ngamma+1):(length(output$phi2))]
      alpha1g <- alpha2g <- vector("list", numBio)
      index = 0
      for(g in 1:numBio){
        alpha1g[[g]] <- alpha1[(index+1):(index + pREvec[g])]
        alpha2g[[g]] <- alpha2[(index+1):(index + pREvec[g])]
        index = index + pREvec[g]
      }
      
      alphaList <- list(alpha1g, alpha2g)
      H01 <- output$H01
      H02 <- output$H02
      
      
      
      if((mvDiffrelative(beta, prebeta, sigma, presigma, gamma1, pregamma1, gamma2, pregamma2,
                 alpha1, prealpha1, alpha2, prealpha2,
                 Sig, preSig, H01, preH01, H02, preH02, tol) == 0) || (iter == maxiter) #|| (!is.list(GetEfun)) || (!is.list(GetMpara))
      ) {
        break
      }
      
    }
    
    if (iter == maxiter) {
      writeLines("program stops because of nonconvergence")
      convergence = 0
      sebeta <- sesigma <- segamma1 <- segamma2 <- sealpha1 <- sealpha2 <- seSig <- vcov <- NULL
      
    } else {
      convergence = 1
      ## run 1 estep here to get posterior mean and covariance --> use estimate to plug in to line47
      
      CUH01 <- rep(0, n)
      CUH02 <- rep(0, n)
      HAZ01 <- rep(0, n)
      HAZ02 <- rep(0, n)
      CumuH01 <- cumsum(H01[, 3])
      CumuH02 <- cumsum(H02[, 3])
      getHazard(CumuH01, CumuH02, survtime, cmprsk, H01, H02, CUH01, CUH02, HAZ01, HAZ02)
      
      data <- list(beta = betaList, gamma1 = gamma1, gamma2 = gamma2,
                   alpha = alphaList, sigma = sigma,
                   Z = subZ, X1 = subX1, Y = subY, Sig = Sig,
                   CUH01 = CUH01,  HAZ01 = HAZ01, CUH02 = CUH02,  HAZ02 = HAZ02,
                   mdataM = mdataM, mdataSM = mdataSM,
                   cmprsk = cmprsk, W = W)
      
      res <- future.apply::future_lapply(seq_len(numSubj), estepMV_worker, future.seed = TRUE, future.scheduling = 2,
                                         data, pREtotal)
      pos.mode <- lapply(res, `[[`, "mode")
      pos.cov  <- lapply(res, function(x) crossprod(x$ccov))
      
      SEest <- getmvCov(beta, gamma1, gamma2,
                          alpha1, alpha2,
                          H01, H02, pos.cov, Sig, sigma, 
                          subX1, subY, subZ, getinit$W, 
                          getinit$survtime,getinit$cmprsk,
                          mdataM, mdataSM, pos.mode)
      
      sebeta <- SEest$sebeta
      sesigma <- SEest$sesigma
      segamma1 <- SEest$segamma1
      segamma2 <- SEest$segamma2
      sealpha1 <- SEest$sealpha1
      sealpha2 <- SEest$sealpha2
      seSig <- SEest$seSig
      vcov <- SEest$vcov
      FisherInfo <- SEest$FisherInfo
      Score <- SEest$Score
      
    }
    
    end_time <- Sys.time()
    
    writeLines("runtime is:")
    print(runtime <- end_time - start_time)
    
    PropComp <- as.data.frame(table(cdata[, survival[2]]))
    call <- match.call()
    
    names(gamma1) <- paste0(namesgamma, "_1")
    names(gamma2) <- paste0(namesgamma, "_2")
    betaList <- vector("list", numBio)
    betacount <- 0
    for (g in 1:numBio) {
      betaList[[g]] <- beta[(betacount+1):(betacount+pbeta[g])]
      names(betaList[[g]]) <- namesbeta[[g]]
      betacount <- betacount + pbeta[g]
    }
    
    names(beta) <- unlist(namesbeta)
    
    result <- list(beta = beta, betaList = output$betaList, gamma1 = gamma1, gamma2 = gamma2, 
                   alpha1 = alpha1, alpha2 = alpha2, H01 = H01, H02 = H02, 
                   Sig = Sig, sigma = sigma, iter = iter, convergence = convergence, 
                   vcov = vcov, FisherInfo = FisherInfo, Score = Score, sebeta = sebeta, segamma1 = segamma1, segamma2 = segamma2, 
                   sealpha1 = sealpha1, sealpha2 = sealpha2, seSig = seSig, sesigma = sesigma, pos.mode = pos.mode, pos.cov = pos.cov,
                   CompetingRisk = CompetingRisk, ydata = rawydata, cdata = rawcdata, 
                   PropEventType = PropComp, LongitudinalSubmodel = long.formula,
                   SurvivalSubmodel = surv.formula, random = random, call = call, id = ID, opt = opt,
                   runtime = runtime)
    
    
    class(result) <- "mvjmcs"
    
    return(result)
    
  } else {
    # ======================== #
    # ~~~~ Single Failure ~~~~ #
    # ======================== #
    
    pREvec <- c()
    
    for(g in 1:numBio){
      pREvec[g] <- ncol(getinit$Z[[g]])
    }
    
    pREtotal <- sum(pREvec)
    
    index = 0
    if (is.null(initial.para)) {
      SigList <- getinit$Sig # need to be 4x4
      Sig <- matrix(rep(0, pREtotal), nrow = pREtotal, ncol = pREtotal)
      for(g in 1:length(pREvec)){
        Sig[(index+1):(index+pREvec[g]),(index+1):(index+pREvec[g])] <- SigList[[g]]
        index = index + pREvec[g]
      }
      
    } else {
      Sig <- getinit$Sig
    }
    
    numSubj <- n
    # opt <- list()
    pos.mode <- vector("list", numSubj)
    subX1 <- subY <- subZ <- vector("list", numSubj)
    pos.cov <- list()
    namesbeta <- vector("list", numBio)
    pbeta <- c()
    for (g in 1:numBio) {
      namesbeta[[g]] <- paste0(names(getinit$beta[[g]]), "_bio", g)
      pbeta[g] <- length(getinit$beta[[g]])
    }
    namesgamma <- names(getinit$gamma1)
    beta <- unlist(getinit$beta)
    betaList <- getinit$beta
    gamma1 <- getinit$gamma1
    alphaList <- getinit$alpha
    alpha1 <- unlist(alphaList[[1]])
    sigma <- getinit$sigma
    
    for(j in 1:numSubj){
      subX1[[j]] <- vector("list", numBio)
      for (g in 1:numBio) {
        
        numRep <- mdataM[[g]][j]
        indexStart <- mdataSM[[g]][j]
        
        subX1[[j]][[g]] <- getinit$X1[[g]][indexStart:(indexStart+numRep-1),, drop = FALSE]
        subY[[j]][[g]] <- getinit$Y[[g]][indexStart:(indexStart+numRep-1)]
        subZ[[j]][[g]] <- getinit$Z[[g]][indexStart:(indexStart+numRep-1),, drop = FALSE]
      }
    }
    
    survtime <- getinit$survtime
    cmprsk <- getinit$cmprsk
    W <- getinit$W
    ngamma <- ncol(W)
    
    old_parallelly <- options(
      parallelly.makeNodePSOCK.rscript_args = c("--vanilla")
    )
    on.exit(options(old_parallelly), add = TRUE)
    old_plan <- future::plan()
    future::plan(future::multisession, workers = cpu.cores)
    on.exit(future::plan(old_plan), add = TRUE)
    
    repeat{
      
      iter <- iter + 1
      
      CUH01 <- rep(0, n) 
      HAZ01 <- rep(0, n) 
      CumuH01 <- cumsum(H01[, 3]) 
      getHazardSF(CumuH01,  survtime, cmprsk, H01, CUH01, HAZ01) 
      
      prebeta <- beta
      pregamma1 <- gamma1
      prealphaList <- alphaList
      prealpha1 <- alpha1 
      preH01 <- H01
      preSig <- Sig
      presigma <- sigma
      
      if (print.para) {
        writeLines("iter is:")
        print(iter)
        writeLines("beta is:")
        print(prebeta)
        writeLines("sigma is:")
        print(presigma)
        writeLines("Sig is:")
        print(preSig)
        writeLines("gamma1 is:")
        print(pregamma1)
        writeLines("alpha1 is:")
        print(prealpha1)
      }
      
      data <- list(beta = betaList, gamma1 = gamma1,
                   alpha = alphaList, sigma = sigma,
                   Z = subZ, X1 = subX1, Y = subY, Sig = Sig,
                   CUH01 = CUH01, HAZ01 = HAZ01, 
                   cmprsk = cmprsk, W = W)
      
      # specify cpu amounts
      
      res <- future.apply::future_lapply(seq_len(numSubj), estepMV_workerSF, future.seed = TRUE, future.scheduling = 2,
                                         data, pREtotal)
      pos.mode <- lapply(res, `[[`, "mode")
      pos.cov <- lapply(res, function(x) crossprod(x$ccov)) 
      
      output <- normalApproxSF(subX1,subY, subZ, W,
                               mdataM, mdataSM,
                               pos.mode,  sigma, pos.cov,
                               H01, survtime, cmprsk,
                               gamma1, alphaList,
                               CUH01,HAZ01, Sig, beta)
      
      beta <- output$beta
      betaList <- output$betaList
      sigma <- output$sigmaVec
      Sig <- output$Sig
      
      gamma1 <- output$phi1[1:ngamma]
      alpha1 <- output$phi1[(ngamma+1):(length(output$phi1))]
      alpha1g <- vector("list", numBio)
      index = 0
      for(g in 1:numBio){
        alpha1g[[g]] <- alpha1[(index+1):(index + pREvec[g])]
        index = index + pREvec[g]
      }
      
      alphaList <- list(alpha1g)
      
      H01 <- output$H01
      
      
      # leave condition
      if((mvDiffrelativeSF(beta, prebeta, sigma, presigma, gamma1, pregamma1, 
                   alpha1, prealpha1,
                   Sig, preSig, H01, preH01, tol) == 0) || (iter == maxiter) #|| (!is.list(GetEfun)) || (!is.list(GetMpara))
      ) {
        break
      }
      
    }
    
    if (iter == maxiter) {
      writeLines("program stops because of nonconvergence")
      convergence = 0
      sebeta <- sesigma <- segamma1 <- sealpha1 <- seSig <- vcov <- NULL
      
    } else {
      convergence = 1
      
      CUH01 <- rep(0, n)
      HAZ01 <- rep(0, n)
      CumuH01 <- cumsum(H01[, 3])
      getHazardSF(CumuH01, survtime, cmprsk, H01, CUH01, HAZ01)
      
      data <- list(beta = betaList, gamma1 = gamma1,
                   alpha = alphaList, sigma = sigma,
                   Z = subZ, X1 = subX1, Y = subY, Sig = Sig,
                   CUH01 = CUH01, HAZ01 = HAZ01,
                   cmprsk = cmprsk, W = W)
      
      res <- future.apply::future_lapply(seq_len(numSubj), estepMV_workerSF, future.seed = TRUE, future.scheduling = 2, data, pREtotal)
      pos.mode <- lapply(res, `[[`, "mode")
      pos.cov  <- lapply(res, function(x) crossprod(x$ccov))
      
      SEest <- getmvCovSF(beta, gamma1,
                          alpha1,  
                          H01, pos.cov, Sig, sigma, 
                          subX1, subY, subZ, getinit$W, 
                          getinit$survtime,getinit$cmprsk,
                          mdataM, mdataSM, pos.mode)
      
      sebeta <- SEest$sebeta
      sesigma <- SEest$sesigma
      segamma1 <- SEest$segamma1
      sealpha1 <- SEest$sealpha1
      seSig <- SEest$seSig
      vcov <- SEest$vcov
      FisherInfo <- SEest$FisherInfo
      Score <- SEest$Score
    }
    
    end_time <- Sys.time()
    print(runtime <- end_time - start_time)
    
    PropComp <- as.data.frame(table(cdata[, survival[2]]))
    call <- match.call()
    
    names(gamma1) <- paste0(namesgamma, "_1")
    betaList <- vector("list", numBio)
    betacount <- 0
    for (g in 1:numBio) {
      betaList[[g]] <- beta[(betacount+1):(betacount+pbeta[g])]
      names(betaList[[g]]) <- namesbeta[[g]]
      betacount <- betacount + pbeta[g]
    }
    
    names(beta) <- unlist(namesbeta)
    
    
    result <- list(beta = beta, betaList = betaList, gamma1 = gamma1, 
                   alpha1 = alpha1, H01 = H01, 
                   Sig = Sig, sigma = sigma, iter = iter, convergence = convergence, 
                   vcov = vcov, FisherInfo = FisherInfo, Score = Score, sebeta = sebeta, segamma1 = segamma1,
                   sealpha1 = sealpha1, seSig = seSig, sesigma = sesigma, 
                   pos.mode = pos.mode, pos.cov = pos.cov,
                   CompetingRisk = CompetingRisk, ydata = rawydata, cdata = rawcdata, 
                   PropEventType = PropComp, LongitudinalSubmodel = long.formula,
                   SurvivalSubmodel = surv.formula, random = random, call = call, id = ID, opt = opt,
                   runtime = runtime)
    
    class(result) <- "mvjmcs"
    
    return(result)
  }
  
}