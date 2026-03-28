GetBayes <- function(beta, sigma, gamma1, gamma2, nu1, nu2, H01, H02, 
                     Sig, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS, initial.optimizer) {
  
  
  n <- nrow(X2)
  p1a <- ncol(Z)
  CUH01 <- rep(0, n)
  CUH02 <- rep(0, n)
  HAZ01 <- rep(0, n)
  HAZ02 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  CumuH02 <- cumsum(H02[, 3])
  
  getHazard(CumuH01, CumuH02, survtime, cmprsk, H01, H02, CUH01, CUH02, HAZ01, HAZ02)
  
  posterior.mode <- matrix(0, nrow = n, ncol = p1a)
  for (i in 1:n) {
    
    if (i != n) {
      subY <- Y[mdataS[i]:(mdataS[i+1]-1)]
      subX1 <- matrix(X1[mdataS[i]:(mdataS[i+1]-1), ], ncol = ncol(X1))
      subZ <- matrix(Z[mdataS[i]:(mdataS[i+1]-1), ], ncol = ncol(Z))
    } else {
      subY <- Y[mdataS[i]:length(Y)]
      subX1 <- matrix(X1[mdataS[i]:length(Y), ], ncol = ncol(X1))
      subZ <- matrix(Z[mdataS[i]:length(Y), ], ncol = ncol(Z))
    }
    CH001 <- CUH01[i]
    CH002 <- CUH02[i]
    HAZ001 <- 1e-3
    HAZ002 <- 1e-3
    
    data <- list(subY, subX1, subZ, t(as.matrix(X2[i, ])), CH001, CH002, 
                 HAZ001, HAZ002, beta, sigma, gamma1, gamma2, nu1, nu2, 
                 Sig, cmprsk[i])
    names(data) <- c("Y", "X", "Z", "X2", "CH01", "CH02", 
                     "HAZ01", "HAZ02", "beta", "sigma",
                     "gamma1", "gamma2", "nu1", "nu2", "Sig", "D")
    opt <- optim(rep(0, p1a), logLikCR.learn, data = data, method = initial.optimizer, hessian = FALSE)
    posterior.mode[i, ] <- opt$par
  }
  
  return(posterior.mode)
  
}

GetBayesSF <- function(beta, sigma, gamma1, nu1, H01, 
                     Sig, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS, initial.optimizer) {
  
  
  n <- nrow(X2)
  p1a <- ncol(Z)
  CUH01 <- rep(0, n)
  HAZ01 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  
  getHazardSF(CumuH01, survtime, cmprsk, H01, CUH01, HAZ01)
  
  posterior.mode <- matrix(0, nrow = n, ncol = p1a)
  for (i in 1:n) {
    
    if (i != n) {
      subY <- Y[mdataS[i]:(mdataS[i+1]-1)]
      subX1 <- matrix(X1[mdataS[i]:(mdataS[i+1]-1), ], ncol = ncol(X1))
      subZ <- matrix(Z[mdataS[i]:(mdataS[i+1]-1), ], ncol = ncol(Z))
    } else {
      subY <- Y[mdataS[i]:length(Y)]
      subX1 <- matrix(X1[mdataS[i]:length(Y), ], ncol = ncol(X1))
      subZ <- matrix(Z[mdataS[i]:length(Y), ], ncol = ncol(Z))
    }
    CH001 <- CUH01[i]
    HAZ001 <- 1e-3
    
    data <- list(subY, subX1, subZ, t(as.matrix(X2[i, ])), CH001,
                 HAZ001, beta, sigma, gamma1, nu1,
                 Sig, cmprsk[i])
    names(data) <- c("Y", "X", "Z", "X2", "CH01", 
                     "HAZ01", "beta", "sigma",
                     "gamma1", "nu1", "Sig", "D")
    opt <- optim(rep(0, p1a), logLik.learn, data = data, method = initial.optimizer, hessian = FALSE)
    posterior.mode[i, ] <- opt$par
  }
  
  return(posterior.mode)
  
}

GetBayesSF.JMH <- function(beta, tau, gamma1, alpha1, vee1, H01, 
                         Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, initial.optimizer) {
  
  
  n <- nrow(X2)
  p1a <- ncol(Z)
  nsig <- p1a + 1
  CUH01 <- rep(0, n)
  HAZ01 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  
  getHazardSF(CumuH01, survtime, cmprsk, H01, CUH01, HAZ01)
  
  posterior.mode <- matrix(0, nrow = n, ncol = nsig)
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
    HAZ001 <- 1e-3
    
    data <- list(subY, subX1, subZ, subW, t(as.matrix(X2[i, ])), CH001, 
                 HAZ001, beta, tau, gamma1, alpha1, vee1, Sig, cmprsk[i])
    names(data) <- c("Y", "X", "Z", "W", "X2", "CH01", "HAZ01", "beta", "tau",
                     "gamma1", "alpha1", "nu1", "Sig", "D")
    opt <- optim(rep(0, nsig), logLik.learn.JMH, data = data, method = initial.optimizer, hessian = FALSE)
    posterior.mode[i, ] <- opt$par
  }
  
  return(posterior.mode)
  
}

GetBayes.JMH <- function(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, 
                     Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, initial.optimizer) {
  
  
  n <- nrow(X2)
  p1a <- ncol(Z)
  nsig <- p1a + 1
  CUH01 <- rep(0, n)
  CUH02 <- rep(0, n)
  HAZ01 <- rep(0, n)
  HAZ02 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  CumuH02 <- cumsum(H02[, 3])
  
  getHazard(CumuH01, CumuH02, survtime, cmprsk, H01, H02, CUH01, CUH02, HAZ01, HAZ02)
  
  posterior.mode <- matrix(0, nrow = n, ncol = nsig)
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
    CH002 <- CUH02[i]
    HAZ001 <- 1e-3
    HAZ002 <- 1e-3
    
    data <- list(subY, subX1, subZ, subW, t(as.matrix(X2[i, ])), CH001, CH002, 
                 HAZ001, HAZ002, beta, tau, gamma1, gamma2, alpha1, alpha2, 
                 vee1, vee2, Sig, cmprsk[i])
    names(data) <- c("Y", "X", "Z", "W", "X2", "CH01", "CH02", 
                     "HAZ01", "HAZ02", "beta", "tau",
                     "gamma1", "gamma2", "alpha1", "alpha2", "nu1", "nu2", "Sig", "D")
    opt <- optim(rep(0, nsig), logLikCR.learn.JMH, data = data, method = initial.optimizer, hessian = FALSE)
    posterior.mode[i, ] <- opt$par
  }
  
  return(posterior.mode)
  
}