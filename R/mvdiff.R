mvDiff <- function(beta, prebeta, sigma, presigma, gamma1, pregamma1, gamma2, pregamma2,
                   alpha1, prealpha1, alpha2, prealpha2,
                   Sig, preSig, H01, preH01, H02, preH02, tol) {
  
  betaAbsdiff <- max(abs(beta - prebeta))
  sigmaAbsdiff <- max(abs(sigma - presigma))
  gamma1Absdiff <- max(abs(gamma1 - pregamma1))
  gamma2Absdiff <- max(abs(gamma2 - pregamma2))
  alpha1Absdiff <- max(abs(alpha1 - prealpha1))
  alpha2Absdiff <- max(abs(alpha2 - prealpha2))
  SigAbsdiff <- max(abs(Sig - preSig))
  H01Absdiff <- max(abs(H01[, 3] - preH01[, 3]))
  H02Absdiff <- max(abs(H02[, 3] - preH02[, 3]))
  
  if ((betaAbsdiff > tol) || (sigmaAbsdiff > tol) || 
      (gamma1Absdiff > tol)|| (gamma2Absdiff > tol) || 
      (alpha1Absdiff > tol) || (alpha2Absdiff > tol) || 
      (SigAbsdiff  > 10*tol) || 
      (H01Absdiff > 10*tol) || (H02Absdiff > 10*tol)) {
    return(1)
  } else {
    return(0)
  }
  
}

mvDiffrelative <- function(beta, prebeta, sigma, presigma, gamma1, pregamma1, gamma2, pregamma2,
                           alpha1, prealpha1, alpha2, prealpha2,
                           Sig, preSig, H01, preH01, H02, preH02, tol) {
  
  betarediff <- max(abs(beta - prebeta)/(abs(prebeta) + 10*tol))
  sigmarediff <- max(abs(sigma - presigma)/(abs(presigma) + 10*tol))
  gamma1rediff <- max(abs(gamma1 - pregamma1)/(abs(pregamma1) + 10*tol))
  gamma2rediff <- max(abs(gamma2 - pregamma2)/(abs(pregamma2) + 10*tol))
  alpha1rediff <- max(abs(alpha1 - prealpha1)/(abs(prealpha1) + 10*tol))
  alpha2rediff <- max(abs(alpha2 - prealpha2)/(abs(prealpha2) + 10*tol))
  Sigrediff <- max(abs(Sig - preSig)/(abs(preSig) + 10*tol))
  H01rediff <- max(abs(H01[, 3] - preH01[, 3])/(abs(preH01[, 3]) + 10*tol))
  H02rediff <- max(abs(H02[, 3] - preH02[, 3])/(abs(preH02[, 3]) + 10*tol))
  
  if ((betarediff > tol) || (sigmarediff > tol) || 
      (gamma1rediff > tol)|| (gamma2rediff > tol) || 
      (alpha1rediff > tol) || (alpha2rediff > tol) || 
      (Sigrediff  > tol) || 
      (H01rediff > tol) || (H02rediff > tol)) {
    return(1)
  } else {
    return(0)
  }
  
}