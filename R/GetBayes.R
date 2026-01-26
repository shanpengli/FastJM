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