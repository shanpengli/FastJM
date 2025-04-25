##' @export
##' 

mvjmcs <- function(ydata, cdata, long.formula,
                         random = NULL, surv.formula,
                         maxiter = 10000, opt = "nlminb", tol = 0.0001, method = c("quad", "noquad"), 
                         model, print.para = TRUE, 
                   initial.para = NULL,
                   quadpoint = 6){
  # "aGH", "normApprox", 
  
  start_time <- Sys.time()
  
  random.form <- all.vars(random)
  RE <- random.form[-length(random.form)]
  ID <- random.form[length(random.form)]
  
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
  
  survival <- all.vars(surv.formula)
  random.form <- all.vars(random)
  ID <- random.form[length(random.form)]
  cnames <- colnames(cdata)
  ynames <- colnames(ydata)
  
  survfmla <- surv.formula
  
  rawydata <- ydata
  rawcdata <- cdata
  
  getinit <- Getmvinit(cdata = cdata, ydata = ydata, long.formula = long.formula,
                       surv.formula = surv.formula,
                       model = "interslope", ID = ID, RE = RE, survinitial = TRUE,
                       REML = TRUE, random = random, opt = "nlminb", initial.para)
  
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
    H01 <- as.matrix(getriskset$tablerisk1)
    
  }
  
  GH.val  <- gauss.quad.prob(quadpoint)
  weight.c <- GH.val$weights # weights
  abscissas.c <- GH.val$nodes # abscissas
  
  survtime <- getinit$survtime
  cmprsk <- getinit$cmprsk
  
  iter=0
  
  
  getriskset <- Getriskset(cdata = getinit$cdata, surv.formula = surv.formula)
  
  ## number of distinct survival time
  H01 <- getriskset$tablerisk1
  H02 <- getriskset$tablerisk2
  
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
      alpha = getinit$alpha,
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
  
  if(method == "quad"){
    output <- getQuadMix(subX1,subY, subZ, getinit$W,
                         mdataM, mdataSM,
                         pos.mode,  getinit$sigma, pos.cov, weight.c, abscissas.c,
                         H01, H02, getinit$survtime, getinit$cmprsk,
                         getinit$gamma1, getinit$gamma2, getinit$alpha,
                         CUH01, CUH02,HAZ01,HAZ02,Sig, subdata$beta)
  }else{
    output <- getNoQuad(subX1,subY, subZ, getinit$W,
                        mdataM, mdataSM,
                        pos.mode,  getinit$sigma, pos.cov, weight.c, abscissas.c,
                        H01, H02, getinit$survtime, getinit$cmprsk,
                        getinit$gamma1, getinit$gamma2, getinit$alpha,
                        CUH01, CUH02,HAZ01,HAZ02,Sig, subdata$beta)
  }
  
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
    
    
    if(method == "quad"){
      output <- getQuadMix(subX1,subY, subZ, getinit$W,
                           mdataM, mdataSM,
                           pos.mode, presigma, pos.cov, weight.c, abscissas.c,
                           H01, H02, getinit$survtime, getinit$cmprsk,
                           data$gamma1, data$gamma2, data$alpha,
                           CUH01, CUH02,HAZ01,HAZ02,preSig, subdata$beta)
      
    }else{
      output <- getNoQuad(subX1,subY, subZ, getinit$W,
                          mdataM, mdataSM,
                          pos.mode, presigma, pos.cov, weight.c, abscissas.c,
                          H01, H02, getinit$survtime, getinit$cmprsk,
                          data$gamma1, data$gamma2, data$alpha,
                          CUH01, CUH02,HAZ01,HAZ02,preSig, subdata$beta)
    }
    
    # PAR UPDATE HERE - for debugging purposes
    beta <- output$beta
    sigma <- output$sigmaVec
    Sig <- output$Sig

    # tempsigma <- rbind(tempsigma, sigma)
    # CHANGE THIS PART
    gamma1 <- output$phi1[1:ngamma]
    gamma2 <- output$phi2[1:ngamma]
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
    SEest = 0
    
  } else{
    SEest <- getmvCov(beta, gamma1, gamma2, 
                    alpha1, alpha2, 
                    H01, H02, pos.cov, Sig, sigma, 
                    subX1, subY, subZ, getinit$W, 
                    getinit$survtime,getinit$cmprsk,
                    mdataM, mdataSM, pos.mode)
  }
  
  end_time <- Sys.time()
  (runtime <- end_time - start_time)
  
  return(list(output = output, re = pos.mode, sigi = pos.cov, 
              beta = beta, sigmaout = sigma, gamma1 = gamma1, gamma2 = gamma2, alpha1 = alpha1, alpha2 = alpha2,
              SEest = SEest, runtime = runtime, iter = iter))
  
}

