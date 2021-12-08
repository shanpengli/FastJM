GetM <- function(GetEfun, beta, gamma1, gamma2, alpha1, alpha2, H01, H02, 
                 Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS) {
  
  if(is.list(GetEfun)) {
    FUNBS <- as.matrix(GetEfun$FUNBS)
    FUNEC <- as.matrix(GetEfun$FUNEC)
    FUNBEC <- as.matrix(GetEfun$FUNBEC)
    FUNBSEC <- as.matrix(GetEfun$FUNBSEC)
    FUNB <- as.matrix(GetEfun$FUNB)
    
    getMpara <- getMC(beta, gamma1, gamma2, alpha1, alpha2, H01, 
                      H02, Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS,
                      FUNBS, FUNEC, FUNBEC, FUNBSEC, FUNB)
      
    return(getMpara)
  
  } else {
    return(0);
  }
}