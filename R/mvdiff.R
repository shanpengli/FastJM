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