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

DiffSFrelative.JMH <- function(beta, prebeta, tau, pretau, gamma1, pregamma1,
                           alpha1, prealpha1, vee1, prevee1,
                           Sig, preSig, H01, preH01, epsilon) {
  
  betarediff <- max(abs(beta - prebeta)/(abs(prebeta) + 10*epsilon))
  taurediff <- max(abs(tau - pretau)/(abs(pretau) + 10*epsilon))
  gamma1rediff <- max(abs(gamma1 - pregamma1)/(abs(pregamma1) + 10*epsilon))
  alpha1rediff <- max(abs(alpha1 - prealpha1)/(abs(prealpha1) + 10*epsilon))
  vee1rediff <- max(abs(vee1 - prevee1)/(abs(prevee1) + 10*epsilon))
  Sigrediff <- max(abs(Sig - preSig)/(abs(preSig) + 10*epsilon))
  H01rediff <- max(abs(H01[, 3] - preH01[, 3])/(abs(preH01[, 3]) + 10*epsilon))
  
  if ((betarediff > epsilon) || (taurediff > epsilon) || (gamma1rediff > epsilon)
      || (alpha1rediff > epsilon) || (vee1rediff > epsilon) || (Sigrediff  > epsilon)
      || (H01rediff > epsilon)) {
    return(1)
  } else {
    return(0)
  }
  
}