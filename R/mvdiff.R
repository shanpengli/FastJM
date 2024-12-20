mvDiff <- function(beta, prebeta, sigma, presigma, gamma1, pregamma1, gamma2, pregamma2,
                   alpha1b1, prealpha1b1, alpha2b1, prealpha2b1,
                   alpha1b2, prealpha1b2, alpha2b2, prealpha2b2,
                   Sig, preSig, H01, preH01, H02, preH02, tol) {
  
  betaAbsdiff <- max(abs(beta - prebeta))
  sigmaAbsdiff <- max(abs(sigma - presigma))
  gamma1Absdiff <- max(abs(gamma1 - pregamma1))
  gamma2Absdiff <- max(abs(gamma2 - pregamma2))
  alpha1b1Absdiff <- max(abs(alpha1b1 - prealpha1b1))
  alpha2b1Absdiff <- max(abs(alpha2b1 - prealpha2b1))
  alpha1b2Absdiff <- max(abs(alpha1b2 - prealpha1b2))
  alpha2b2Absdiff <- max(abs(alpha2b2 - prealpha2b2))
  SigAbsdiff <- max(abs(Sig - preSig))
  H01Absdiff <- max(abs(H01[, 3] - preH01[, 3]))
  H02Absdiff <- max(abs(H02[, 3] - preH02[, 3]))
  
  if ((betaAbsdiff > tol) || (sigmaAbsdiff > tol) || (gamma1Absdiff > tol)
      || (gamma2Absdiff > tol) || (alpha1b1Absdiff > tol) || (alpha2b1Absdiff > tol)
      || (alpha1b2Absdiff > tol) || (alpha2b2Absdiff > tol)
      || (SigAbsdiff  > 10*tol) || (H01Absdiff > 10*tol) || (H02Absdiff > 10*tol)) {
    return(1)
  } else {
    return(0)
  }
  
}