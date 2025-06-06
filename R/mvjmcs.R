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
##'
##' @return A list containing:
##' \item{output}{EM algorithm output from final iteration}
##' \item{re}{Estimated random effects for each subject}
##' \item{sigi}{Estimated random effects covariance matrices (posterior) for each subject}
##' \item{beta}{Estimated fixed effects for longitudinal models}
##' \item{sigmaout}{Biomarker error variance estimates}
##' \item{gamma1}{Fixed effects for cause 1}
##' \item{gamma2}{Fixed effects for cause 2 (if \code{CompetingRisk = TRUE})}
##' \item{alpha1}{Association parameters for cause 1}
##' \item{alpha2}{Association parameters for cause 2}
##' \item{SEest}{Estimated standard errors of all parameters}
##' \item{runtime}{Elapsed run time}
##' \item{iter}{Number of EM iterations run}
##'
##' @examples
##' 
##' 
##'   require(FastJM)
##'   require(survival
##' 
##'   # Fit joint model with two biomarkers
##'   fit <-mvjmcs(ydata, cdata, long.formula = list(Y1 ~ X11 + X12 + time, Y2 ~ X11 + X12 + time),
##'                 random = list(~time| ID, ~1|ID),
##'                 surv.formula =Surv(survtime, cmprsk) ~ X21 + X22, maxiter = 50, opt = "nlminb", tol = 0.001)
##'   fit
##'   
##'   # Obtain the variance-variance matrix of all parameter estimates
##'   



##' @seealso \code{\link{jmcs}}, \code{\link{survfitjmcs}}, \code{\link{AUCjmcs}}, \code{\link{MAEQjmcs}}, \code{\link{PEjmcs}}
##' @export

mvjmcs <- function(ydata, cdata, long.formula,
                   random = NULL, surv.formula,
                   maxiter = 10000, opt = "nlminb", tol = 0.005,
                   print.para = TRUE, 
                   initial.para = NULL, Gamma1Cal, Gamma2Cal){
  
  start_time <- Sys.time()
  
  #
  # if(numBio == 1){
  #   mvjmcsSB()
  # }
  

  
  if(is.list(long.formula)){
    numBio = length(long.formula)
  } else {
    numBio = 1
  }
  
  random.form <- RE  <- model <- vector("list", numBio)
  
  # <- ID
  
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
      # ID[[g]] <- random.form[[g]][length(all.vars(random[[g]]))]
    }
    ID <- random.form[[g]][length(all.vars(random[[g]]))]
  } else {
    for(g in 1:numBio){
      random.form[[g]] <- all.vars(random)
      if(length(random) == 1){
        # RE[[g]] <- NULL
        model[[g]] <- "intercept"
      } else {
        RE[[g]] <- random.form[[g]][-length(all.vars(random))]
        model[[g]] <- "interslope"
      }
      # ID[[g]] <- random.form[[g]][length(all.vars(random))]
    }
    ID <- random.form[[g]][length(all.vars(random))]
  }
  
  lengthb <- length(long.formula)
  long <- list()
  
  # check for every biomarker
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
  
  # survival <- all.vars(surv.formula)
  # random.form <- all.vars(random)
  # ID <- random.form[length(random.form)]
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
  
  # need to rearrange this part
  numBio = length(long.formula)
  
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
  
  if(CompetingRisk == TRUE){
    
  
  
  CUH01 <- rep(0, n)
  CUH02 <- rep(0, n)
  HAZ01 <- rep(0, n)
  HAZ02 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  CumuH02 <- cumsum(H02[, 3])
  
  getHazard(CumuH01, CumuH02, getinit$survtime, getinit$cmprsk, H01, H02, CUH01, CUH02, HAZ01, HAZ02)
  
  HAZ0 <- list(HAZ01, HAZ02)
  
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
  
  data <- list(beta = getinit$beta, gamma1 = getinit$gamma1, gamma2 = getinit$gamma2,
               alpha = getinit$alpha, sigma = getinit$sigma,
               Z = getinit$Z, X1 = getinit$X1, Y = getinit$Y, Sig = Sig,
               CUH01 = CUH01, CUH02 = CUH02, HAZ01 = HAZ01, HAZ02 = HAZ02,
               mdataM = mdataM, mdataSM = mdataSM,
               cmprsk = getinit$cmprsk, W = getinit$W)
  
  numSubj <- n
  numBio <- length(data$X1)
  opt <- list()
  pos.mode <- vector("list", numSubj)
  subX1 <- subY <- subZ <- vector("list", numSubj)
  pos.cov <- list()
  submdataM <- vector("list", numBio)
  submdataSM <- vector("list", numBio)
  
  for(j in 1:numSubj) {
    subX1[[j]] <- vector("list", numBio)
    for (g in 1:numBio) {
      
      numRep <- data$mdataM[[g]][j]
      indexStart <- data$mdataSM[[g]][j]
      
      subX1[[j]][[g]] <- data$X1[[g]][indexStart:(indexStart+numRep-1),, drop = FALSE]
      subY[[j]][[g]] <- data$Y[[g]][indexStart:(indexStart+numRep-1)]
      subZ[[j]][[g]] <- data$Z[[g]][indexStart:(indexStart+numRep-1),, drop = FALSE]
    }
    
    subCUH01 <- data$CUH01[j]
    subCUH02 <- data$CUH02[j]
    subHAZ01 <- data$HAZ01[j]
    subHAZ02 <- data$HAZ02[j]
    subcmprsk <- data$cmprsk[j]
    subW <- t(data$W[j, ])
    
    subdata <- list(
      beta = data$beta,
      gamma1 = data$gamma1,
      gamma2 = data$gamma2,
      alpha = data$alpha,
      sigma = data$sigma,
      Z = subZ[[j]],
      X1 = subX1[[j]],
      Y = subY[[j]],
      Sig = data$Sig,
      CUH01 = subCUH01,
      CUH02 = subCUH02,
      HAZ01 = subHAZ01,
      HAZ02 = subHAZ02,
      mdataM = submdataM,
      mdataSM = submdataSM,
      cmprsk = subcmprsk,
      W = subW
    )
    
    opt <- optim(
      par = c(rep(0, pREtotal)), # CHANGE THIS PART
      getbSig,
      getbSig_grad,
      data = subdata,
      method = "BFGS",
      hessian = TRUE
    )
    pos.mode[[j]] <- opt$par
    pos.cov[[j]] <- solve(opt$hessian)
    
  }
  
  survtime <- getinit$survtime
  cmprsk <- getinit$cmprsk
  
  # if(method == "quad"){
  #   output <- getQuadMix(subX1,subY, subZ, getinit$W,
  #                        mdataM, mdataSM,
  #                        pos.mode,  getinit$sigma, pos.cov, weight.c, abscissas.c,
  #                        H01, H02, getinit$survtime, getinit$cmprsk,
  #                        getinit$gamma1, getinit$gamma2, getinit$alpha,
  #                        CUH01, CUH02,HAZ01,HAZ02,Sig, subdata$beta)
  # } else {
  output <- normalApprox(subX1,subY, subZ, getinit$W,
                         mdataM, mdataSM,
                         pos.mode,  getinit$sigma, pos.cov,
                         H01, H02, getinit$survtime, getinit$cmprsk,
                         getinit$gamma1, getinit$gamma2, getinit$alpha,
                         CUH01, CUH02,HAZ01,HAZ02, Sig, subdata$beta)
  # }
  
  # PAR UPDATE HERE
  
  tempbeta <- output$beta
  # tempsigma <- c(output$sigma1, output$sigma2)
  tempphi1 <- output$phi1
  tempphi2 <- output$phi2
  
  it <- 1
  tempSig <- list()
  tempSig[[it]] <- output$Sig
  ngamma <- ncol(data$W)
  index = 0
  gamma1 <- output$phi1[(index+1):ngamma]
  gamma2 <- output$phi2[(index+1):ngamma]
  
  # gamma1 hiim<- c(1,0.5)
  # gamma2 <- c(-0.5,0.5)
  index = index + ngamma
  #NEED TO ADJUST THIS PART
  alpha1 <- output$phi1[(ngamma+1):(length(output$phi1))] # 3 and 2 come from competing risk
  alpha2 <- output$phi2[(ngamma+1):(length(output$phi2))]
  
  alpha1g <- alpha2g <- vector("list", numBio)
  index = 0
  for(g in 1:numBio){
    alpha1g[[g]] <- alpha1[(index+1):(index + pREvec[g])]
    alpha2g[[g]] <- alpha2[(index+1):(index + pREvec[g])]
    index = index + pREvec[g]
  }
  
  
  alphaList <- list(alpha1g, alpha2g)
  
  iter=0
  beta <- output$beta
  sigma <- output$sigmaVec
  Sig <- output$Sig

  repeat{
    
    iter <- iter + 1
    prebeta <- beta
    pregamma1 <- gamma1
    pregamma2 <- gamma2
    
    prealphaList <- alphaList
    prealpha1 <- alpha1 # 3 and 2 come from competing risk
    prealpha2 <- alpha2
    preH01 <- H01
    preH02 <- H02
    preSig <- Sig
    presigma <- sigma # for checking
    
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
    
    CUH01 <- rep(0, n)
    CUH02 <- rep(0, n)
    HAZ01 <- rep(0, n)
    HAZ02 <- rep(0, n)
    
    CumuH01 <- cumsum(output$H01[, 3])
    CumuH02 <- cumsum(output$H02[, 3])
    
    getHazard(CumuH01, CumuH02, survtime, cmprsk, preH01, preH02, CUH01, CUH02, HAZ01, HAZ02)
    
    data <- list(beta = prebeta, gamma1 = pregamma1, gamma2 = pregamma2,
                 alpha = prealphaList, sigma = presigma,
                 Z = getinit$Z, X1 = getinit$X1, Y = getinit$Y, Sig = preSig,
                 CUH01 = CUH01, CUH02 = CUH02, HAZ01 = HAZ01, HAZ02 = HAZ02,
                 mdataM = mdataM, mdataSM = mdataSM,
                 cmprsk = getinit$cmprsk, W = getinit$W)
    
    for(j in 1:numSubj) {
      
      subCUH01 <- CUH01[j]
      subCUH02 <- CUH02[j]
      subHAZ01 <- HAZ01[j]
      subHAZ02 <- HAZ02[j]
      subcmprsk <- data$cmprsk[j]
      subW <- t(data$W[j, ])
      
      subdata <- list(
        beta = output$betaList,
        gamma1 = data$gamma1,
        gamma2 = data$gamma2,
        alpha = data$alpha,
        sigma = data$sigma,
        Z = subZ[[j]],
        X1 = subX1[[j]],
        Y = subY[[j]],
        Sig = data$Sig,
        CUH01 = subCUH01,
        CUH02 = subCUH02,
        HAZ01 = subHAZ01,
        HAZ02 = subHAZ02,
        mdataM = submdataM,
        mdataSM = submdataSM,
        cmprsk = subcmprsk,
        W = subW
      )
      
      opt <- optim(
        par = rep(0, pREtotal),
        getbSig,
        getbSig_grad,
        data = subdata,
        method = "BFGS",
        hessian = TRUE
      )
      
      pos.mode[[j]] <- opt$par
      pos.cov[[j]] <- solve(opt$hessian)
    }
    
    
    # if(method == "quad"){
    #   output <- getQuadMix(subX1,subY, subZ, getinit$W,
    #                        mdataM, mdataSM,
    #                        pos.mode, presigma, pos.cov, weight.c, abscissas.c,
    #                        H01, H02, getinit$survtime, getinit$cmprsk,
    #                        data$gamma1, data$gamma2, data$alpha,
    #                        CUH01, CUH02,HAZ01,HAZ02,preSig, subdata$beta)
    #   
    # } else {
    output <- normalApprox(subX1,subY, subZ, getinit$W,
                           mdataM, mdataSM,
                           pos.mode, presigma, pos.cov,
                           H01, H02, getinit$survtime, getinit$cmprsk,
                           data$gamma1, data$gamma2, data$alpha,
                           CUH01, CUH02,HAZ01,HAZ02,preSig, subdata$beta)
    # }
    
    # PAR UPDATE HERE - for debugging purposes
    beta <- output$beta
    sigma <- output$sigmaVec
    Sig <- output$Sig
    
    # tempsigma <- rbind(tempsigma, sigma)
    # CHANGE THIS PART
    gamma1 <- output$phi1[1:ngamma]
    gamma2 <- output$phi2[1:ngamma]
    # gamma1 <- c(1,0.5)
    # gamma2 <- c(-0.5,0.5)
    alpha1 <- output$phi1[(ngamma+1):(length(output$phi1))] # 3 and 2 come from competing risk
    alpha2 <- output$phi2[(ngamma+1):(length(output$phi2))]
    alpha1g <- alpha2g <- vector("list", numBio)
    index = 0
    for(g in 1:numBio){
      alpha1g[[g]] <- alpha1[(index+1):(index + pREvec[g])]
      alpha2g[[g]] <- alpha2[(index+1):(index + pREvec[g])]
      index = index + pREvec[g]
    }
    
    alphaList <- list(alpha1g, alpha2g)
    
    it <- it + 1
    tempSig[[it]] <- Sig
    H01 <- output$H01
    H02 <- output$H02
    
    # leave condition
    if((mvDiff(beta, prebeta, sigma, presigma, gamma1, pregamma1, gamma2, pregamma2,
               alpha1, prealpha1, alpha2, prealpha2,
               Sig, preSig, H01, preH01, H02, preH02, tol) == 0) || (iter == maxiter) #|| (!is.list(GetEfun)) || (!is.list(GetMpara))
    ) {
      break
    }
    
  }
  
  if (iter == maxiter) {
    writeLines("program stops because of nonconvergence")
    convergence = 0
    sebeta <- sesigma <- segamma1 <- segamma2 <- sealpha1 <- sealpha2 <- seSig <- NULL
    
  } else {
    convergence = 1
    SEest <- getmvCov(beta, gamma1, gamma2,
                      alpha1, alpha2,
                      H01, H02, pos.cov, Sig, sigma,
                      subX1, subY, subZ, getinit$W,
                      getinit$survtime,getinit$cmprsk,
                      mdataM, mdataSM, pos.mode)
    
    # SEest <- getmvCovTest(beta, gamma1, gamma2,
    #                   alpha1, alpha2,
    #                   H01, H02, pos.cov, Sig, sigma,
    #                   subX1, subY, subZ, getinit$W,
    #                   getinit$survtime,getinit$cmprsk,
    #                   mdataM, mdataSM, pos.mode, Gamma1Cal, Gamma2Cal)
    
    
    
    sebeta <- SEest$sebeta
    sesigma <- SEest$sesigma
    segamma1 <- SEest$segamma1
    segamma2 <- SEest$segamma2
    sealpha1 <- SEest$sealpha1
    sealpha2 <- SEest$sealpha2
    seSig <- SEest$seSig
    vcov <- SEest$vcov
    
    
  }
  
  end_time <- Sys.time()
  print(runtime <- end_time - start_time)
  
  # return(list(output = output, re = pos.mode, sigi = pos.cov, 
  #             beta = beta, sigmaout = sigma, gamma1 = gamma1, gamma2 = gamma2, alpha1 = alpha1, alpha2 = alpha2,
  #             SEest = SEest, runtime = runtime, iter = iter))
  
  
  return(list(beta = beta, gamma1 = gamma1, gamma2 = gamma2, 
              alpha1 = alpha1, alpha2 = alpha2, H01 = H01, H02 = H02, 
              Sig = Sig, sigma = sigma, iter = iter, convergence = convergence, 
              vcov = vcov, sebeta = sebeta, segamma1 = segamma1, segamma2 = segamma2, 
              sealpha1 = sealpha1, sealpha1 = sealpha2, seSig = seSig, sesigma = sesigma, pos.mode = pos.mode, 
              pos.cov = pos.cov, runtime = runtime))
  } else {
  # ~~~~~~~~~~~~~
  # SINGLE FAILURE
  # ~~~~~~~~~~~~~
    
    CUH01 <- rep(0, n)
    HAZ01 <- rep(0, n)
    CumuH01 <- cumsum(H01[, 3])
    
    getHazardSF(CumuH01, getinit$survtime, getinit$cmprsk, H01,  CUH01, HAZ01)
    
    HAZ0 <- list(HAZ01)
    
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
    
    data <- list(beta = getinit$beta, gamma1 = getinit$gamma1,
                 alpha = getinit$alpha, sigma = getinit$sigma,
                 Z = getinit$Z, X1 = getinit$X1, Y = getinit$Y, Sig = Sig,
                 CUH01 = CUH01, HAZ01 = HAZ01, 
                 mdataM = mdataM, mdataSM = mdataSM,
                 cmprsk = getinit$cmprsk, W = getinit$W)
    
    numSubj <- n
    numBio <- length(data$X1)
    opt <- list()
    pos.mode <- vector("list", numSubj)
    subX1 <- subY <- subZ <- vector("list", numSubj)
    pos.cov <- list()
    submdataM <- vector("list", numBio)
    submdataSM <- vector("list", numBio)
    
    for(j in 1:numSubj) {
      subX1[[j]] <- vector("list", numBio)
      for (g in 1:numBio) {
        
        numRep <- data$mdataM[[g]][j]
        indexStart <- data$mdataSM[[g]][j]
        
        subX1[[j]][[g]] <- data$X1[[g]][indexStart:(indexStart+numRep-1),, drop = FALSE]
        subY[[j]][[g]] <- data$Y[[g]][indexStart:(indexStart+numRep-1)]
        subZ[[j]][[g]] <- data$Z[[g]][indexStart:(indexStart+numRep-1),, drop = FALSE]
      }
      
      subCUH01 <- data$CUH01[j]
      subHAZ01 <- data$HAZ01[j]
      subcmprsk <- data$cmprsk[j]
      subW <- t(data$W[j, ])
      
      subdata <- list(
        beta = data$beta,
        gamma1 = data$gamma1,
        alpha = getinit$alpha,
        sigma = data$sigma,
        Z = subZ[[j]],
        X1 = subX1[[j]],
        Y = subY[[j]],
        Sig = data$Sig,
        CUH01 = subCUH01,
        HAZ01 = subHAZ01,
        mdataM = submdataM,
        mdataSM = submdataSM,
        cmprsk = subcmprsk,
        W = subW
      )
      
      opt <- optim(
        par = c(rep(0, pREtotal)), # CHANGE THIS PART
        getbSigSF,
        getbSig_gradSF,
        data = subdata,
        method = "BFGS",
        hessian = TRUE
      )
      pos.mode[[j]] <- opt$par
      pos.cov[[j]] <- solve(opt$hessian)
      
    }
    
    survtime <- getinit$survtime
    cmprsk <- getinit$cmprsk
    
    # if(method == "quad"){
    #   output <- getQuadMix(subX1,subY, subZ, getinit$W,
    #                        mdataM, mdataSM,
    #                        pos.mode,  getinit$sigma, pos.cov, weight.c, abscissas.c,
    #                        H01, H02, getinit$survtime, getinit$cmprsk,
    #                        getinit$gamma1, getinit$gamma2, getinit$alpha,
    #                        CUH01, CUH02,HAZ01,HAZ02,Sig, subdata$beta)
    # } else {
    output <- normalApproxSF(subX1,subY, subZ, getinit$W,
                           mdataM, mdataSM,
                           pos.mode,  getinit$sigma, pos.cov, 
                           H01, getinit$survtime, getinit$cmprsk,
                           getinit$gamma1, getinit$alpha,
                           CUH01,HAZ01,Sig, subdata$beta)
    # }
    
    # PAR UPDATE HERE
    
    tempbeta <- output$beta
    # tempsigma <- c(output$sigma1, output$sigma2)
    tempphi1 <- output$phi1
    
    it <- 1
    tempSig <- list()
    tempSig[[it]] <- output$Sig
    ngamma <- ncol(data$W)
    index = 0
    gamma1 <- output$phi1[(index+1):ngamma]
    
    # gamma1 hiim<- c(1,0.5)
    # gamma2 <- c(-0.5,0.5)
    index = index + ngamma
    #NEED TO ADJUST THIS PART
    alpha1 <- output$phi1[(ngamma+1):(length(output$phi1))] # 3 and 2 come from competing risk
    
    alpha1g <- vector("list", numBio)
    index = 0
    for(g in 1:numBio){
      alpha1g[[g]] <- alpha1[(index+1):(index + pREvec[g])]
      index = index + pREvec[g]
    }
    
    
    alphaList <- list(alpha1g)
    
    iter=0
    beta <- output$beta
    sigma <- output$sigmaVec
    Sig <- output$Sig
    
    repeat{
      
      iter <- iter + 1
      prebeta <- beta
      pregamma1 <- gamma1
      prealphaList <- alphaList
      prealpha1 <- alpha1
      preH01 <- H01
      preSig <- Sig
      presigma <- sigma # for checking
      
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
      
      CUH01 <- rep(0, n)
      HAZ01 <- rep(0, n)
      
      CumuH01 <- cumsum(output$H01[, 3])
      
      getHazardSF(CumuH01, survtime, cmprsk, preH01, CUH01, HAZ01)
      
      data <- list(beta = prebeta, gamma1 = pregamma1,
                   alpha = prealphaList, sigma = presigma,
                   Z = getinit$Z, X1 = getinit$X1, Y = getinit$Y, Sig = preSig,
                   CUH01 = CUH01, HAZ01 = HAZ01,
                   mdataM = mdataM, mdataSM = mdataSM,
                   cmprsk = getinit$cmprsk, W = getinit$W)
      
      for(j in 1:numSubj) {
        
        subCUH01 <- CUH01[j]
        subHAZ01 <- HAZ01[j]
        subcmprsk <- data$cmprsk[j]
        subW <- t(data$W[j, ])
        
        subdata <- list(
          beta = output$betaList,
          gamma1 = data$gamma1,
          alpha = data$alpha,
          sigma = data$sigma,
          Z = subZ[[j]],
          X1 = subX1[[j]],
          Y = subY[[j]],
          Sig = data$Sig,
          CUH01 = subCUH01,
          HAZ01 = subHAZ01,
          mdataM = submdataM,
          mdataSM = submdataSM,
          cmprsk = subcmprsk,
          W = subW
        )
        
        opt <- optim(
          par = rep(0, pREtotal),
          getbSigSF,
          getbSig_gradSF,
          data = subdata,
          method = "BFGS",
          hessian = TRUE
        )
        
        pos.mode[[j]] <- opt$par
        pos.cov[[j]] <- solve(opt$hessian)
      }
      
      output <- normalApproxSF(subX1,subY, subZ, getinit$W,
                             mdataM, mdataSM,
                             pos.mode, presigma, pos.cov,
                             H01, getinit$survtime, getinit$cmprsk,
                             data$gamma1, data$alpha,
                             CUH01,HAZ01, preSig, subdata$beta)
      # }
      
      # PAR UPDATE HERE - for debugging purposes
      beta <- output$beta
      sigma <- output$sigmaVec
      Sig <- output$Sig
      
      # tempsigma <- rbind(tempsigma, sigma)
      # CHANGE THIS PART
      gamma1 <- output$phi1[1:ngamma]
      alpha1 <- output$phi1[(ngamma+1):(length(output$phi1))] # 3 and 2 come from competing risk
      alpha1g <- vector("list", numBio)
      index = 0
      for(g in 1:numBio){
        alpha1g[[g]] <- alpha1[(index+1):(index + pREvec[g])]
        index = index + pREvec[g]
      }
      
      alphaList <- list(alpha1g)
      
      it <- it + 1
      tempSig[[it]] <- Sig
      H01 <- output$H01
      
      # leave condition
      if((mvDiffSF(beta, prebeta, sigma, presigma, gamma1, pregamma1, 
                 alpha1, prealpha1,
                 Sig, preSig, H01, preH01, tol) == 0) || (iter == maxiter) #|| (!is.list(GetEfun)) || (!is.list(GetMpara))
      ) {
        break
      }
      
    }
    
    if (iter == maxiter) {
      writeLines("program stops because of nonconvergence")
      convergence = 0
      sebeta <- sesigma <- segamma1 <- sealpha1 <- seSig <- NULL
      
    } else {
      convergence = 1
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
    }
    
    end_time <- Sys.time()
    print(runtime <- end_time - start_time)
    
    
    return(list(beta = beta, gamma1 = gamma1, 
                alpha1 = alpha1, H01 = H01, 
                Sig = Sig, sigma = sigma, iter = iter, convergence = convergence, 
                vcov = vcov, sebeta = sebeta, segamma1 = segamma1,
                sealpha1 = sealpha1, seSig = seSig, sesigma = sesigma, pos.mode = pos.mode, pos.cov = pos.cov,
                runtime = runtime))
  }
  
}
