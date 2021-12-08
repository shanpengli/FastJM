DiffSF <- function(beta, prebeta, sigma, presigma, gamma1, pregamma1,
                 alpha1, prealpha1, Sig, preSig, H01, preH01, epsilon) {
  
  betaAbsdiff <- max(abs(beta - prebeta))
  sigmaAbsdiff <- max(abs(sigma - presigma))
  gamma1Absdiff <- max(abs(gamma1 - pregamma1))
  alpha1Absdiff <- max(abs(alpha1 - prealpha1))
  SigAbsdiff <- max(abs(Sig - preSig))
  H01Absdiff <- max(abs(H01[, 3] - preH01[, 3]))
  
  if ((betaAbsdiff > epsilon) || (sigmaAbsdiff > epsilon) || (gamma1Absdiff > epsilon)
       || (alpha1Absdiff > epsilon) || (SigAbsdiff  > epsilon) || (H01Absdiff > 10*epsilon)) {
    return(1)
  } else {
    return(0)
  }
  
}