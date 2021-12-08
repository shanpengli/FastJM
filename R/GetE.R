GetE <- function(beta, gamma1, gamma2, alpha1, alpha2, H01, H02, 
                 Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix,
                 method, Posbi, Poscov) {
  
  
  n <- nrow(X2)
  p1a <- ncol(Z)
  
  CUH01 <- rep(0, n)
  CUH02 <- rep(0, n)
  HAZ01 <- rep(0, n)
  HAZ02 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  CumuH02 <- cumsum(H02[, 3])
  
  getHazard(CumuH01, CumuH02, survtime, cmprsk, H01, H02, CUH01, CUH02, HAZ01, HAZ02)
  
  if (method == "standard") {
    status = getECstandard(beta, gamma1, gamma2, alpha1, alpha2,
                           Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, 
                   wsmatrix, CUH01, CUH02, HAZ01, HAZ02)
  } else {
    status = getECpseudo(beta, gamma1, gamma2, alpha1, alpha2,
                           Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix,
                           wsmatrix, CUH01, CUH02, HAZ01, HAZ02, Posbi, Poscov)

  }
  
  
  return(status)
  
}