GetESF <- function(beta, gamma1, alpha1, H01, Sig, sigma, Z, X1, Y, X2, survtime, 
                   cmprsk, mdata, mdataS, xsmatrix, wsmatrix, method, Posbi, Poscov) {
  
  
  n <- nrow(X2)
  p1a <- ncol(Z)
  
  CUH01 <- rep(0, n)
  HAZ01 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  
  getHazardSF(CumuH01, survtime, cmprsk, H01, CUH01, HAZ01)
  
  if (method == "standard") {
    status = getECstandardSF(beta, gamma1, alpha1, Sig, sigma, 
                             Z, X1, Y, X2, survtime, cmprsk, mdata, 
                             mdataS, xsmatrix, wsmatrix, CUH01, HAZ01)
  } else {
    status = getECpseudoSF(beta, gamma1, alpha1, Sig, sigma, 
                             Z, X1, Y, X2, survtime, cmprsk, mdata, 
                             mdataS, xsmatrix, wsmatrix, CUH01, HAZ01, Posbi, Poscov)
  }
  
  return(status)
  
}

GetESF.JMH <- function(beta, tau, gamma1, alpha1, vee1, H01, 
                   Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix) {
  
  
  n <- nrow(X2)
  p1a <- ncol(Z)
  
  CUH01 <- rep(0, n)
  HAZ01 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  
  getHazardSF(CumuH01, survtime, cmprsk, H01, CUH01, HAZ01)
  
  
  status = getECSF(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, 
                   survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, HAZ01)
  
  return(status)
  
}

GetESFad.JMH <- function(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, 
                     survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, initial.optimizer) {
  
  
  n <- nrow(X2)
  p1a <- ncol(Z)
  nsig <- p1a + 1
  CUH01 <- rep(0, n)
  HAZ01 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  
  getHazardSF(CumuH01, survtime, cmprsk, H01, CUH01, HAZ01)
  
  posterior.mode <- matrix(0, nrow = n, ncol = nsig)
  posterior.var <- matrix(0, nrow = nsig*n, ncol = nsig)
  for (i in 1:n) {
    
    if (i != n) {
      subY <- Y[mdataS[i]:(mdataS[i+1]-1)]
      subW <- matrix(W[mdataS[i]:(mdataS[i+1]-1), ], ncol = ncol(W))
      subX1 <- matrix(X1[mdataS[i]:(mdataS[i+1]-1), ], ncol = ncol(X1))
      subZ <- matrix(Z[mdataS[i]:(mdataS[i+1]-1), ], ncol = ncol(Z))
    } else {
      subY <- Y[mdataS[i]:length(Y)]
      subW <- matrix(W[mdataS[i]:length(Y), ], ncol = ncol(W))
      subX1 <- matrix(X1[mdataS[i]:length(Y), ], ncol = ncol(X1))
      subZ <- matrix(Z[mdataS[i]:length(Y), ], ncol = ncol(Z))
    }
    CH001 <- CUH01[i]
    HAZ001 <- HAZ01[i]
    
    data <- list(subY, subX1, subZ, subW, t(as.matrix(X2[i, ])), CH001, 
                 HAZ001, beta, tau, gamma1, alpha1, 
                 vee1, Sig, cmprsk[i])
    names(data) <- c("Y", "X", "Z", "W", "X2", "CH01", 
                     "HAZ01", "beta", "tau",
                     "gamma1", "alpha1", "nu1", "Sig", "D")
    opt <- optim(rep(0, nsig), logLik.learn.JMH, data = data, method = initial.optimizer, hessian = TRUE)
    posterior.mode[i, ] <- opt$par
    posterior.var[(nsig*(i-1) + 1):(i * nsig), 1:nsig] <- solve(opt$hessian)
  }
  
  status = getECSFad(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, 
                     survtime, cmprsk, mdata, mdataS, xsmatrix,
                     wsmatrix, CUH01, HAZ01, posterior.mode, posterior.var)
  
  return(status)
  
}