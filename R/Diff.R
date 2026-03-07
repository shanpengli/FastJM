Diff <- function(beta, prebeta, sigma, presigma, gamma1, pregamma1, gamma2, pregamma2,
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
  
  if ((betaAbsdiff > tol) || (sigmaAbsdiff > tol) || (gamma1Absdiff > tol)
      || (gamma2Absdiff > tol) || (alpha1Absdiff > tol) || (alpha2Absdiff > tol)
      || (SigAbsdiff  > 10*tol) || (H01Absdiff > 10*tol) || (H02Absdiff > 10*tol)) {
    return(1)
  } else {
    return(0)
  }
  
}

Diffrelative.JMH <- function(beta, prebeta, tau, pretau, gamma1, pregamma1, gamma2, pregamma2,
                         alpha1, prealpha1, alpha2, prealpha2, vee1, prevee1, vee2, prevee2,
                         Sig, preSig, H01, preH01, H02, preH02, epsilon) {
  
  betarediff <- max(abs(beta - prebeta)/(abs(prebeta) + 10*epsilon))
  taurediff <- max(abs(tau - pretau)/(abs(pretau) + 10*epsilon))
  gamma1rediff <- max(abs(gamma1 - pregamma1)/(abs(pregamma1) + 10*epsilon))
  gamma2rediff <- max(abs(gamma2 - pregamma2)/(abs(pregamma2) + 10*epsilon))
  alpha1rediff <- max(abs(alpha1 - prealpha1)/(abs(prealpha1) + 10*epsilon))
  alpha2rediff <- max(abs(alpha2 - prealpha2)/(abs(prealpha2) + 10*epsilon))
  vee1rediff <- max(abs(vee1 - prevee1)/(abs(prevee1) + 10*epsilon))
  vee2rediff <- max(abs(vee2 - prevee2)/(abs(prevee2) + 10*epsilon))
  Sigrediff <- max(abs(Sig - preSig)/(abs(preSig) + 10*epsilon))
  H01rediff <- max(abs(H01[, 3] - preH01[, 3])/(abs(preH01[, 3]) + 10*epsilon))
  H02rediff <- max(abs(H02[, 3] - preH02[, 3])/(abs(preH02[, 3]) + 10*epsilon))
  
  if ((betarediff > epsilon) || (taurediff > epsilon) || (gamma1rediff > epsilon)
      || (gamma2rediff > epsilon) || (alpha1rediff > epsilon) || (alpha2rediff > epsilon)
      || (vee1rediff > epsilon) || (vee2rediff > epsilon) || (Sigrediff  > epsilon)
      || (H01rediff > epsilon) || (H02rediff > epsilon)) {
    return(1)
  } else {
    return(0)
  }
  
}