mvDiffSF <- function(beta, prebeta, sigma, presigma, gamma1, pregamma1,
                   alpha1, prealpha1,
                   Sig, preSig, H01, preH01, tol) {
  
  betaAbsdiff <- max(abs(beta - prebeta))
  sigmaAbsdiff <- max(abs(sigma - presigma))
  gamma1Absdiff <- max(abs(gamma1 - pregamma1))
  alpha1Absdiff <- max(abs(alpha1 - prealpha1))
  SigAbsdiff <- max(abs(Sig - preSig))
  H01Absdiff <- max(abs(H01[, 3] - preH01[, 3]))
  
  if ((betaAbsdiff > tol) || (sigmaAbsdiff > tol) || 
      (gamma1Absdiff > tol)|| 
      (alpha1Absdiff > tol) || 
      (SigAbsdiff  > 10*tol) || 
      (H01Absdiff > 10*tol) ) {
    return(1)
  } else {
    return(0)
  }
  
}

mvDiffrelativeSF <- function(beta, prebeta, sigma, presigma, gamma1, pregamma1,
                           alpha1, prealpha1, Sig, preSig, H01, preH01, tol) {
  
  betarediff <- max(abs(beta - prebeta)/(abs(prebeta) + 10*tol))
  sigmarediff <- max(abs(sigma - presigma)/(abs(presigma) + 10*tol))
  gamma1rediff <- max(abs(gamma1 - pregamma1)/(abs(pregamma1) + 10*tol))
  alpha1rediff <- max(abs(alpha1 - prealpha1)/(abs(prealpha1) + 10*tol))
  Sigrediff <- max(abs(Sig - preSig)/(abs(preSig) + 10*tol))
  H01rediff <- max(abs(H01[, 3] - preH01[, 3])/(abs(preH01[, 3]) + 10*tol))
  
  if ((betarediff > tol) || (sigmarediff > tol) || 
      (gamma1rediff > tol) || (alpha1rediff > tol) || 
      (Sigrediff  > tol) || (H01rediff > tol)) {
    return(1)
  } else {
    return(0)
  }
  
}