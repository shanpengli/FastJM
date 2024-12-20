##' @export
##' 

runEMbeta <- function(ydata, cdata, long.formula, 
                      random = NULL, surv.formula, 
                      maxiter = 10000, opt = "nlminb", tol = 0.0001, method = c("aGH", "normApprox"), model,
                      print.para = TRUE){
  
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
                     REML = TRUE, random = random, opt = "nlminb")
  
  ydata <- getinit$ydata # need to rearrange this part
  cdata <- getinit$cdata
  
  mdataS1 <- mdataS2 <- c()
  mdata1 <- getinit$mdata[[1]]
  n1 <- nrow(mdata1)
  mdata1 <- as.data.frame(mdata1)
  mdata1 <- as.vector(mdata1$ni)
  mdataS1[1] <- 1
  mdataCum1 <- cumsum(mdata1)
  mdata21 <- mdata1 - 1
  mdataS1[2:n1] <- mdataCum1[2:n1] - mdata21[2:n1]
  
  mdata2 <- getinit$mdata[[2]]
  n2 <- nrow(mdata2)
  mdata2 <- as.data.frame(mdata2)
  mdata2 <- as.vector(mdata2$ni)
  mdataS2[1] <- 1
  mdataCum2 <- cumsum(mdata2)
  mdata22 <- mdata1 - 1
  mdataS2[2:n2] <- mdataCum2[2:n2] - mdata21[2:n2]
  
  mdataM <- list(mdata1, mdata2)
  mdataSM <- list(mdataS1, mdataS2)
  
  
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
  
  # get weights - could make this into a function
  GH.val  <- gauss.quad.prob(10)
  weight.c <- GH.val$weights # ws matrix
  abscissas.c <- GH.val$nodes # xs matrix
  
  quadpoint = 6
  gq_vals <- statmod::gauss.quad(n = 6, kind = "hermite")
  xs <- gq_vals$nodes
  ws <- gq_vals$weights
  
  
  xsmatrix <- matrix(0, nrow = 4, ncol = quadpoint^4)
  wsmatrix <- xsmatrix
  
  xsmatrix[4, ] <- rep(xs, quadpoint^3)
  
  
  
  Total <- NULL
  for (i in 1:quadpoint) {
    sub <- rep(xs[i], quadpoint)
    Total <- c(Total, sub)
  }
  xsmatrix[3, ] <- rep(Total, quadpoint^2)
  
  Total <- NULL
  for (i in 1:quadpoint) {
    sub <- rep(xs[i], quadpoint^2)
    Total <- c(Total, sub)
  }
  xsmatrix[2, ] <- rep(Total, quadpoint)
  
  Total <- NULL
  for (i in 1:quadpoint) {
    sub <- rep(xs[i], quadpoint^3)
    Total <- c(Total, sub)
  }
  xsmatrix[1, ] <- Total
  
  xsmatrix <- t(xsmatrix)
  
  
  wsmatrix[4, ] <- rep(ws, quadpoint^3)
  
  Total <- NULL
  for (i in 1:quadpoint) {
    sub <- rep(ws[i], quadpoint)
    Total <- c(Total, sub)
  }
  wsmatrix[3, ] <- rep(Total, quadpoint^2)
  
  Total <- NULL
  for (i in 1:quadpoint) {
    sub <- rep(ws[i], quadpoint^2)
    Total <- c(Total, sub)
  }
  wsmatrix[2, ] <- rep(Total, quadpoint)
  Total <- NULL
  for (i in 1:quadpoint) {
    sub <- rep(ws[i], quadpoint^3)
    Total <- c(Total, sub)
  }
  wsmatrix[1, ] <- Total
  wsmatrix <- t(wsmatrix)
  
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
  
  SigList <- getinit$Sig # need to be 4x4
  Sig11 <- SigList[[1]]
  Sig22<- SigList[[2]] # first biomarker/ #second biomarker
  HAZ0 <- list(HAZ01, HAZ02)
  
  p1a <- ncol(getinit$Z[[1]]) # dim of random effect
  p2a <- ncol(getinit$Z[[2]])
  # Sig <- matrix(rep(0, p1a + p2a), nrow = p1a + p2a, ncol = p1a + p2a)
  # Sig[1:p1a,1:p1a] <- Sig11
  # Sig[(p1a+1):(p1a + p2a),(p1a+1):(p1a+p2a)] <- Sig22
  # 
  # 
  # data <- list(beta = getinit$beta, gamma1 = getinit$gamma1, gamma2 = getinit$gamma2, 
  #              alpha = getinit$alpha, sigma = getinit$sigma,
  #              Z = getinit$Z, X1 = getinit$X1, Y = getinit$Y, Sig = Sig,
  #              CUH01 = CUH01, CUH02 = CUH02, HAZ01 = HAZ01, HAZ02 = HAZ02,
  #              mdataM = mdataM, mdataSM = mdataSM, 
  #              cmprsk = getinit$cmprsk, W = getinit$W)
  
  
  beta = list(beta1 = c(5, 1.5, 2, 1), beta2 = c(10, 1, 2, 1))
  sigma = list(0.5, 0.5)
  gamma1 = c(1, 0.5)
  gamma2 = c(-0.5, 0.5)
  Sig <- diag(1,4)
  alpha = list(alpha1 = list(alpha11 = c(0.5, 0.7),
                             alpha12 = c(-0.5, 0.5)),
               alpha2 = list(alpha21 = c(0.5, 0.7),
                             alpha22 = c(-0.5, 0.5)))
  
  data <- list(beta = getinit$beta, gamma1 = gamma1, gamma2 = gamma2, 
               alpha = alpha, sigma = sigma,
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
  
  
  # GH.val  <- statmod::gauss.quad.prob(10)
  # weight.c <- GH.val$weights # ws matrix
  # abscissas.c <- GH.val$nodes # xs matrix
  
  # quadpoint = 6
  # gq_vals <- statmod::gauss.quad(n = 6, kind = "hermite")
  # xs <- gq_vals$nodes
  # ws <- gq_vals$weights
  # 
  # 
  # xsmatrix <- matrix(0, nrow = 4, ncol = quadpoint^4)
  # wsmatrix <- xsmatrix
  # 
  # xsmatrix[4, ] <- rep(xs, quadpoint^3)
  # 
  # 
  # 
  # Total <- NULL
  # for (i in 1:quadpoint) {
  #   sub <- rep(xs[i], quadpoint)
  #   Total <- c(Total, sub)
  # }
  # xsmatrix[3, ] <- rep(Total, quadpoint^2)
  # 
  # Total <- NULL
  # for (i in 1:quadpoint) {
  #   sub <- rep(xs[i], quadpoint^2)
  #   Total <- c(Total, sub)
  # }
  # xsmatrix[2, ] <- rep(Total, quadpoint)
  # 
  # Total <- NULL
  # for (i in 1:quadpoint) {
  #   sub <- rep(xs[i], quadpoint^3)
  #   Total <- c(Total, sub)
  # }
  # xsmatrix[1, ] <- Total
  # 
  # xsmatrix <- t(xsmatrix)
  # 
  # 
  # wsmatrix[4, ] <- rep(ws, quadpoint^3)
  # 
  # Total <- NULL
  # for (i in 1:quadpoint) {
  #   sub <- rep(ws[i], quadpoint)
  #   Total <- c(Total, sub)
  # }
  # wsmatrix[3, ] <- rep(Total, quadpoint^2)
  # 
  # Total <- NULL
  # for (i in 1:quadpoint) {
  #   sub <- rep(ws[i], quadpoint^2)
  #   Total <- c(Total, sub)
  # }
  # wsmatrix[2, ] <- rep(Total, quadpoint)
  # Total <- NULL
  # for (i in 1:quadpoint) {
  #   sub <- rep(ws[i], quadpoint^3)
  #   Total <- c(Total, sub)
  # }
  # wsmatrix[1, ] <- Total
  # wsmatrix <- t(wsmatrix)
  
  survtime <- getinit$survtime
  cmprsk <- getinit$cmprsk
  # output <- testC(subX1,subY, subZ, data$W,
  #                 mdataM, mdataSM, 
  #                 pos.mode, getinit$sigma, pos.var, weight.c, abscissas.c, 
  #                 H01, H02, survtime,cmprsk,
  #                 getinit$gamma1, getinit$gamma2, getinit$alpha)
  
  # output <- testC(subX1,subY, subZ, data$W,
  #                 mdataM, mdataSM, 
  #                 pos.mode, getinit$sigma, pos.var, weight.c, abscissas.c, 
  #                 H01, H02,getinit$survtime, getinit$cmprsk, 
  #                 getinit$gamma1, getinit$gamma2, getinit$alpha, xsmatrix, wsmatrix,
  #                 CUH01, CUH02,HAZ01,HAZ02,Sig, 10, getinit$beta)  
  
  
  ## temp keeping everything else constant
  
  
  
  
  # output <- testC(subX1,subY, subZ, data$W,
  #                 mdataM, mdataSM,
  #                 pos.mode, sigma, pos.var, weight.c, abscissas.c,
  #                 H01, H02,getinit$survtime, getinit$cmprsk,
  #                 gamma1, gamma2, alpha, xsmatrix, wsmatrix,
  #                 CUH01, CUH02,HAZ01,HAZ02,Sig, getinit$beta)
  
  if(method == "aGH"){
    output <- getMaGH(subX1,subY, subZ, getinit$W,
                      mdataM, mdataSM, 
                      pos.mode, sigma, pos.var, weight.c, abscissas.c, 
                      H01, H02, getinit$survtime, getinit$cmprsk, 
                      gamma1, gamma2, alpha, xsmatrix, wsmatrix,
                      CUH01, CUH02,HAZ01,HAZ02,Sig, subdata$beta)
  }else{
    output <- getMNA(subX1,subY, subZ, getinit$W,
                     mdataM, mdataSM, 
                     pos.mode, sigma, pos.var, weight.c, abscissas.c, 
                     H01, H02, getinit$survtime, getinit$cmprsk, 
                     gamma1, gamma2, alpha, xsmatrix, wsmatrix,
                     CUH01, CUH02,HAZ01,HAZ02,Sig, subdata$beta)
  }
  
  
  tempbeta <- output$beta
  
  n <- nrow(data$W)
  p1a <- ncol(data$Z[[1]])
  p2a <- ncol(data$Z[[1]])
  
  # gamma1 <- output$phi1[1:2]
  # gamma2 <- output$phi2[1:2]
  alphaList <- list(list(output$phi1[3:4], output$phi2[3:4]),list(output$phi1[5:6], output$phi2[5:6]))
  sigmaList <- list(output$sigma1, output$sigma2)
  #sigmaList <- list(output$sigma1q, output$sigma2q)
  iter=0
  
  
  start_time <- Sys.time()
  
  repeat{
    iter <- iter + 1
    prebeta <- output$beta
    
    if (print.para) {
      writeLines("iter is:")
      print(iter)
      writeLines("beta is:")
      print(output$beta)
      # writeLines("gamma1 is:")
      # print(gamma1)
      # writeLines("gamma2 is:")
      # print(gamma2)
      # writeLines("nu1 is:")
      # print(alpha1)
      # writeLines("nu2 is:")
      # print(alpha2)
      # writeLines("Sig is:")
      # print(Sig)
      # writeLines("Error variance is:")
      # print(sigma)
    }
    
    # pregamma1 <- gamma1
    # pregamma2 <- gamma2
    # prealphaList <- alphaList
    # prealpha1b1 <- prealphaList[[1]][[1]]
    # prealpha1b2 <- prealphaList[[1]][[2]]
    # prealpha2b1 <- prealphaList[[2]][[1]]
    # prealpha2b2 <- prealphaList[[2]][[2]]
    # preH01 <- output$H01
    # preH02 <- output$H02
    # preSig <- output$Sig
    # presigmaList <-  c(output$sigma1, output$sigma2)
    # #presigmaList <- c(output$sigma1q, output$sigma2q)
    # preH01 <- output$H01
    # preH02 <- output$H02
    # #preH01 <- output$H01q
    # #preH02 <- output$H02q
    
    
    # CUH01 <- rep(0, n)
    # CUH02 <- rep(0, n)
    # HAZ01 <- rep(0, n)
    # HAZ02 <- rep(0, n)
    # 
    # CumuH01 <- cumsum(output$H01[, 3])
    # CumuH02 <- cumsum(output$H02[, 3])
    # 
    # getHazard(CumuH01, CumuH02, survtime, cmprsk, preH01, preH02, CUH01, CUH02, HAZ01, HAZ02)
    
    
    for(j in 1:numSubj) {
      
      subCUH01 <- data$CUH01[j]
      subCUH02 <- data$CUH02[j]
      subHAZ01 <- data$HAZ01[j]
      subHAZ02 <- data$HAZ02[j]
      subcmprsk <- data$cmprsk[j]
      subW <- data$W[j, ]
      
      # subdata <- list(
      #   beta = list(output$beta1, output$beta2),
      #   gamma1 = pregamma1,
      #   gamma2 = pregamma2,
      #   alpha = prealphaList,
      #   sigma = presigmaList,
      #   Z = subZ[[j]],
      #   X1 = subX1[[j]],
      #   Y = subY[[j]],
      #   Sig = preSig,
      #   CUH01 = subCUH01,
      #   CUH02 = subCUH02,
      #   HAZ01 = subHAZ01,
      #   HAZ02 = subHAZ02,
      #   mdataM = submdataM,
      #   mdataSM = submdataSM,
      #   cmprsk = subcmprsk,
      #   W = subW
      # )
      
      subdata <- list(
        #beta = output$beta,
        beta = list(output$beta1, output$beta2),
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
    
    gamma1 = c(1, 0.5)
    gamma2 = c(-0.5, 0.5)
    
    
    
    
    # output <- testC(subX1,subY, subZ, data$W,
    #       mdataM, mdataSM, 
    #       pos.mode, presigmaList, pos.var, weight.c, abscissas.c, 
    #       preH01, preH02, survtime,cmprsk,
    #       pregamma1, pregamma2, prealphaList)
    
    # output <- testC(subX1,subY, subZ, getinit$W,
    #                 mdataM, mdataSM, 
    #                 pos.mode, getinit$sigma, pos.var, weight.c, abscissas.c, 
    #                 H01, H02,getinit$survtime, getinit$cmprsk, 
    #                 getinit$gamma1, getinit$gamma2, getinit$alpha, xsmatrix, wsmatrix,
    #                 CUH01, CUH02,HAZ01,HAZ02,Sig,10, subdata$beta)  
    
    # temp keeping everything else constant
    # output <- testC(subX1,subY, subZ, getinit$W,
    #                 mdataM, mdataSM, 
    #                 pos.mode, sigma, pos.var, weight.c, abscissas.c, 
    #                 H01, H02, getinit$survtime, getinit$cmprsk, 
    #                 gamma1, gamma2, alpha, xsmatrix, wsmatrix,
    #                 CUH01, CUH02,HAZ01,HAZ02,Sig, subdata$beta)  
    
    
    if(method == "aGH"){
      output <- getMaGH(subX1,subY, subZ, getinit$W,
                        mdataM, mdataSM, 
                        pos.mode, sigma, pos.var, weight.c, abscissas.c, 
                        H01, H02, getinit$survtime, getinit$cmprsk, 
                        gamma1, gamma2, alpha, xsmatrix, wsmatrix,
                        CUH01, CUH02,HAZ01,HAZ02,Sig, subdata$beta)
    }else{
      output <- getMNA(subX1,subY, subZ, getinit$W,
                       mdataM, mdataSM, 
                       pos.mode, sigma, pos.var, weight.c, abscissas.c, 
                       H01, H02, getinit$survtime, getinit$cmprsk, 
                       gamma1, gamma2, alpha, xsmatrix, wsmatrix,
                       CUH01, CUH02,HAZ01,HAZ02,Sig, subdata$beta)
    }
    
    tempbeta <- cbind(tempbeta, output$beta)
    
    beta <- output$beta
    #gamma1 <- output$phi1[1:2]
    #gamma2 <- output$phi2[1:2]
    #alphaList <- list(list(output$phi1[3:4], output$phi2[3:4]),list(output$phi1[5:6], output$phi2[5:6]))
    #alpha1b1 <- alphaList[[1]][[1]]
    # alpha2b1 <- alphaList[[1]][[2]]
    # alpha1b2 <- alphaList[[2]][[1]]
    # alpha2b2 <- alphaList[[2]][[2]]
    # sigma <- c(output$sigma1, output$sigma2)
    # sigma <- c(output$sigma1q, output$sigma2q)
    #Sig <- output$Sig
    #H01 <- output$H01
    #H02 <- output$H02
    #H01 <- output$H01q
    #H02 <- output$H02q
    
    
    #sigmaList <- list(output$sigma1, output$sigma2)
    sigmaList <- list(output$sigma1q, output$sigma2q)
    
    # if((Diff(beta, prebeta, sigma, presigmaList, gamma1, pregamma1, gamma2, pregamma2,
    #          alpha1b1, prealpha1b1, alpha2b1, prealpha2b1, 
    #          alpha1b2, prealpha1b2, alpha2b2, prealpha2b2,
    #          Sig, preSig, H01, preH01, H02, preH02, tol) == 0) || (iter == maxiter) #|| (!is.list(GetEfun)) || (!is.list(GetMpara))
    #    ) {
    
    if(max(abs(beta-prebeta) < tol) || (iter == maxiter)){
      break
    }
    
  }
  
  if (iter == maxiter) {
    writeLines("program stops because of nonconvergence")
    convergence = 0
  }
  
  end_time <- Sys.time()
  (runtime <- end_time - start_time)
  return(list(output = output, re = pos.mode, sigi = pos.var, betaout = tempbeta, runtime = runtime))
  
}

