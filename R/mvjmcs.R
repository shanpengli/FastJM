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
  mdataM <- mdataSM <- vector("list", numBio)
  
  for(g in 1:numBio){
    ydata <- getinit$ydata[[g]]
    mdata <- getinit$mdata[[g]]
    n <- nrow(mdata)
    mdata <- as.data.frame(mdata)
    mdata <- as.vector(mdata$ni)
    mdataS[1] <- 1
    mdataCum <- cumsum(mdata)
    mdata2 <- mdata - 1
    mdataS[2:n] <- mdataCum[2:n] - mdata2[2:n]
  }
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
  
  CUH01 <- rep(0, n1)
  CUH02 <- rep(0, n1)
  HAZ01 <- rep(0, n1)
  HAZ02 <- rep(0, n1)
  
  CumuH01 <- cumsum(H01[, 3])
  CumuH02 <- cumsum(H02[, 3])
  
  getHazard(CumuH01, CumuH02, getinit$survtime, getinit$cmprsk, H01, H02, CUH01, CUH02, HAZ01, HAZ02)
  
  HAZ0 <- list(HAZ01, HAZ02)
  
  p1a <- ncol(getinit$Z[[1]]) # dim of random effect
  p2a <- ncol(getinit$Z[[2]])
  
  if (is.null(initial.para)) {
    
    SigList <- getinit$Sig # need to be 4x4
    Sig11 <- SigList[[1]]
    Sig22<- SigList[[2]] # first biomarker/ #second biomarker
    Sig <- matrix(rep(0, p1a + p2a), nrow = p1a + p2a, ncol = p1a + p2a)
    Sig[1:p1a,1:p1a] <- Sig11
    Sig[(p1a+1):(p1a + p2a),(p1a+1):(p1a+p2a)] <- Sig22
    
  } else {
    Sig <- getinit$Sig
  }

  
  
  # -----------------------
  # PARAMETERS DEFINED HERE - FOR DEBUGGING
  # -----------------------
  
  #beta = list(beta1 = c(5, 1.5, 2, 1), beta2 = c(10, 1, 2, 1))
  #sigma = list(0.5, 0.5)
  #gamma1 = c(1, 0.5)
  #gamma2 = c(-0.5, 0.5)
  #Sig <- diag(1,4)
  # alpha = list(alpha1 = list(alpha11 = c(0.5, 0.7),
  #                            alpha12 = c(-0.5, 0.5)),
  #              alpha2 = list(alpha21 = c(0.5, 0.7),
  #                            alpha22 = c(-0.5, 0.5)))
  
  
  data <- list(beta = getinit$beta, gamma1 = getinit$gamma1, gamma2 = getinit$gamma2,
               alpha = getinit$alpha, sigma = getinit$sigma,
               Z = getinit$Z, X1 = getinit$X1, Y = getinit$Y, Sig = Sig,
               CUH01 = CUH01, CUH02 = CUH02, HAZ01 = HAZ01, HAZ02 = HAZ02,
               mdataM = mdataM, mdataSM = mdataSM,
               cmprsk = getinit$cmprsk, W = getinit$W)
  
  
  numSubj <- length(mdata1)
  numBio <- length(data$X1)
  opt <- list()
  pos.mode <- vector("list", numSubj)
  subX1 <- subY <- subZ <- vector("list", numSubj)
  pos.var <- list()
  submdataM <- vector("list", numBio)
  submdataSM <- vector("list", numBio)
  
  for(j in 1:numSubj) {
    for (g in 1:numBio) {
      numSubj <- length(data$mdataM[[g]])
      
      numRep <- data$mdataM[[g]][j]
      indexStart <- data$mdataSM[[g]][j]
      
      subX1[[j]][[g]] <- matrix(nrow = numRep, ncol = ncol(data$X1[[g]]))
      subY[[j]][[g]] <- matrix(nrow = numRep, ncol = 1)
      subZ[[j]][[g]] <- matrix(nrow = numRep, ncol = ncol(data$Z[[g]]))
      
      
      # change to vector
      submdataM[g] <- data$mdataM[[g]][j]
      submdataSM[g] <- data$mdataSM[[g]][j]
      for (i in 1:numRep) {
        subX1[[j]][[g]][i, ] <-  data$X1[[g]][indexStart + i - 1, ]
        subY[[j]][[g]][i, ] <- data$Y[[g]][indexStart + i - 1]
        subZ[[j]][[g]][i, ] <- data$Z[[g]][indexStart + i - 1, ]
      }
      
    }
    
    subCUH01 <- data$CUH01[j]
    subCUH02 <- data$CUH02[j]
    subHAZ01 <- data$HAZ01[j]
    subHAZ02 <- data$HAZ02[j]
    subcmprsk <- data$cmprsk[j]
    subW <- data$W[j, ]
    
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
      par = c(rep(0, 4)),
      getbSig,
      data = subdata,
      method = "BFGS",
      hessian = TRUE
    )
    
    # pos.mode[[j]] <- matrix(nrow = p1a + p2a, ncol = 1)
    # pos.var[[j]] <- matrix(nrow = p1a + p2a, ncol = p1a + p2a)
    
    pos.mode[[j]][[1]] <- opt$par[1:p1a]
    pos.mode[[j]][[2]] <- opt$par[(p1a+1):(p1a+p2a)]
    pos.var[[j]] <- solve(opt$hessian)
    
    
    
  }
  
  
  
  
  survtime <- getinit$survtime
  cmprsk <- getinit$cmprsk
  
  if(method == "quad"){
    output <- getQuadMix(subX1,subY, subZ, getinit$W,
                         mdataM, mdataSM,
                         pos.mode, getinit$sigma, pos.var, weight.c, abscissas.c,
                         H01, H02, getinit$survtime, getinit$cmprsk,
                         getinit$gamma1, getinit$gamma2, getinit$alpha,
                         CUH01, CUH02,HAZ01,HAZ02,Sig, subdata$beta)
  }else{
    output <- getNoQuad(subX1,subY, subZ, getinit$W,
                        mdataM, mdataSM,
                        pos.mode, getinit$sigma, pos.var, weight.c, abscissas.c,
                        H01, H02, getinit$survtime, getinit$cmprsk,
                        getinit$gamma1, getinit$gamma2, getinit$alpha,
                        CUH01, CUH02,HAZ01,HAZ02,Sig, subdata$beta)
  }
  
  # PAR UPDATE HERE
  
  tempbeta <- output$beta
  tempsigma <- c(output$sigma1, output$sigma2)
  tempphi1 <- output$phi1
  tempphi2 <- output$phi2
  
  it <- 1
  tempSig <- list()
  tempSig[[it]] <- output$Sig
  
  
  
  n <- nrow(data$W)
  p1a <- ncol(data$Z[[1]])
  p2a <- ncol(data$Z[[1]])
  
  gamma1 <- output$phi1[1:2]
  gamma2 <- output$phi2[1:2]
  alphaList <- list(list(output$phi1[3:4], output$phi1[5:6]),list(output$phi2[3:4], output$phi2[5:6]))
  sigmaList <- list(output$sigma1, output$sigma2)
  
  iter=0
  beta <- output$beta
  

  
  repeat{
    
    iter <- iter + 1
    
    
    # if(iter == 2){
    #   writeLines("n is:")
    #   print(j)
    # }
    
    
    prebeta <- beta
    
    pregamma1 <- gamma1
    pregamma2 <- gamma2
    prealphaList <- alphaList
    prealpha1b1 <- prealphaList[[1]][[1]]
    prealpha1b2 <- prealphaList[[1]][[2]]
    prealpha2b1 <- prealphaList[[2]][[1]]
    prealpha2b2 <- prealphaList[[2]][[2]]
    
    # prealpha1b1 <- prealphaList[[1]][(1:2)]
    # prealpha1b2 <- prealphaList[[1]][(3:4)]
    # prealpha2b1 <- prealphaList[[2]][(1:2)]
    # prealpha2b2 <- prealphaList[[2]][(3:4)]
    preH01 <- output$H01
    preH02 <- output$H02
    preSig <- output$Sig
    presigmaList <-  list(output$sigma1, output$sigma2)
    presigma <- c(output$sigma1, output$sigma2) # for checking
    
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
      print(c(prealpha1b1, prealpha1b2))
      writeLines("alpha2 is:")
      print(c(prealpha2b1, prealpha2b2))
      # writeLines("Sig is:")
      # print(Sig)
      # writeLines("Error variance is:")
      # print(sigma)
    }
    
    
    CUH01 <- rep(0, n)
    CUH02 <- rep(0, n)
    HAZ01 <- rep(0, n)
    HAZ02 <- rep(0, n)
    
    CumuH01 <- cumsum(output$H01[, 3])
    CumuH02 <- cumsum(output$H02[, 3])
    
    getHazard(CumuH01, CumuH02, survtime, cmprsk, preH01, preH02, CUH01, CUH02, HAZ01, HAZ02)
    
    data <- list(beta = prebeta, gamma1 = pregamma1, gamma2 = pregamma2,
                 alpha = prealphaList, sigma = presigmaList,
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
      subW <- data$W[j, ]
      
      
      subdata <- list(
        #beta = output$beta,
        beta = list(output$beta1, output$beta2),
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
        par = c(pos.mode[[j]][[1]], pos.mode[[j]][[2]]),
        getbSig,
        data = subdata,
        method = "BFGS",
        hessian = TRUE
      )
      
      
      pos.mode[[j]][[1]] <- opt$par[1:p1a]
      pos.mode[[j]][[2]] <- opt$par[(p1a+1):(p1a+p2a)]
      pos.var[[j]] <- solve(opt$hessian)
      
      
    }
    
    
    if(method == "quad"){
      output <- getQuadMix(subX1,subY, subZ, getinit$W,
                           mdataM, mdataSM,
                           pos.mode, presigmaList, pos.var, weight.c, abscissas.c,
                           H01, H02, getinit$survtime, getinit$cmprsk,
                           data$gamma1, data$gamma2, data$alpha,
                           CUH01, CUH02,HAZ01,HAZ02,preSig, subdata$beta)
      
    }else{
      output <- getNoQuad(subX1,subY, subZ, getinit$W,
                          mdataM, mdataSM,
                          pos.mode, presigmaList, pos.var, weight.c, abscissas.c,
                          H01, H02, getinit$survtime, getinit$cmprsk,
                          data$gamma1, data$gamma2, data$alpha,
                          CUH01, CUH02,HAZ01,HAZ02,preSig, subdata$beta)
    }
    
    # PAR UPDATE HERE - for debugging purposes
    
    tempbeta <- cbind(tempbeta, output$beta)
    beta <- output$beta
    sigma <- c(output$sigma1, output$sigma2)
    Sig <- output$Sig
    tempphi1 <- c(tempphi1, output$phi1)
    tempphi2 <- c(tempphi2, output$phi2)
    tempsigma <- rbind(tempsigma, sigma)
    gamma1 <- output$phi1[1:2]
    gamma2 <- output$phi2[1:2]
    alphaList <- list(list(output$phi1[3:4], output$phi1[5:6]),list(output$phi2[3:4], output$phi2[5:6]))
    alpha1b1 <- alphaList[[1]][[1]]
    alpha1b2 <- alphaList[[1]][[2]]
    alpha2b1 <- alphaList[[2]][[1]]
    alpha2b2 <- alphaList[[2]][[2]]
    it <- it + 1
    tempSig[[it]] <- Sig
    H01 <- output$H01
    H02 <- output$H02
    
    
    sigmaList <- list(output$sigma1, output$sigma2)
    
    
    # leave condition
    if((mvDiff(beta, prebeta, sigma, presigma, gamma1, pregamma1, gamma2, pregamma2,
             alpha1b1, prealpha1b1, alpha2b1, prealpha2b1,
             alpha1b2, prealpha1b2, alpha2b2, prealpha2b2,
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
                    c(alpha1b1, alpha1b2), c(alpha2b1, alpha2b2), 
                    H01, H02, pos.var, Sig, sigma, 
                    subX1, subY, subZ, getinit$W, 
                    getinit$survtime,getinit$cmprsk,
                    mdata, mdataSM, pos.mode)
  }
  
  end_time <- Sys.time()
  (runtime <- end_time - start_time)
  
  return(list(output = output, re = pos.mode, sigi = pos.var, 
              betaout = tempbeta, sigmaout = tempsigma, 
              phi1out = tempphi1, phi2out = tempphi2, SEest,
              runtime = runtime))
  
}