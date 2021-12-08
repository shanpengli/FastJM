getLoglikeSF <- function(beta, gamma1, alpha1, H01, 
                       Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, 
                       mdataS, xsmatrix, wsmatrix, method, Posbi, Pscov) {
  
  
  n <- nrow(X2)
  p1a <- ncol(Z)
  
  CUH01 <- rep(0, n)
  HAZ01 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  
  getHazardSF(CumuH01, survtime, cmprsk, H01, CUH01, HAZ01)
  
  if (method == "standard") {
    status = getloglikeCstandardSF(beta, gamma1, alpha1,
                                 Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, 
                                 mdata, mdataS, xsmatrix, wsmatrix, 
                                 CUH01, HAZ01)
  } else {
    status = getloglikeCpseudoSF(beta, gamma1, alpha1,
                                   Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, 
                                   mdata, mdataS, xsmatrix, wsmatrix, 
                                   CUH01, HAZ01, Posbi, Pscov)
  }
  
  
  return(status)
  
}