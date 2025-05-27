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