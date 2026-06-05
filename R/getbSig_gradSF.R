getbSig_gradSF <- function(bSig, data){
  
  
  # sigma <- exp(data$W %*% data$tau + w)  # calculates sigma; don't need assume homogeneous error variance
  # sum over G
  # Sig: covariance matrix of B
  # sigma: vector of error variance for each biomarker
  # unique to each value
  
  latAsso <- data$latAsso
  Xs_i <- data$Xs_i
  Zs_i <- data$Zs_i
  Y <- data$Y
  X <- data$X # update so both biomarkers accountted for
  Z <- data$Z
  
  beta <- data$beta
  
  alphaList <- data$alphaList
  sigma <- data$sigma
  
  Sig <- data$Sig
  
  # Sig11 <- SigList[[1]]
  # Sig22<- SigList[[2]] # first biomarker/ #second biomarker
  
  # survival
  CH01 <- data$CH01
  HAZ01 <- data$HAZ01
  Wcmprsk <- data$Wcmprsk
  Wx <- as.matrix(data$Wx)
  
  if(ncol(Wx) ==1){
    Wx <- t(Wx)
  }
  
  gamma1 <- as.matrix(data$gamma1) # vector
  
  # p12 <- nrow(Z[[1]])
  # p22<- nrow(Z[[2]])
  
  
  # turn into list/loop
  pREvec <- c()
  if(is.list(Z)){
    for(g in 1:length(Z)){
      pREvec[g] <- ncol(Z[[g]])
    }
  } else {
    pREvec[1] <- ncol(Z)
  }
  
  q = sum(pREvec)
  
  index = 0
  b <- vector("list", length(pREvec))
  bfull <- c()
  for(p in 1:length(pREvec)){
    b[[p]] <- matrix(bSig[(index + 1):(index + pREvec[p])], nrow = pREvec[p], ncol = 1)
    bfull[(index + 1):(index + pREvec[p])] <- t(b[[p]])
    index = index + pREvec[p]
  }
  
  bfull <- as.matrix(bfull)
  
  total <- c()
  sum.alpha1i <- 0
  
  # need to generalize here
  dEta1 <- numeric(q)
  
  # longitudinal portion
  index <- 0
  for (g in 1:length(Y)) {
    
    Yi <- as.matrix(Y[[g]])
    Xi <- as.matrix(X[[g]])
    
    if(is.list(beta)){
      betai <- as.matrix(beta[[g]])
    }else{
      betai <- as.matrix(beta)
    }
    
    Zi <- as.matrix(Z[[g]])
    bi <- as.matrix(b[[g]])
    sigmai <- sigma[g]
    alpha1 <- alphaList[[1]] # risk 1
    # gets for each biomarker
    
    if(is.list(alpha1)){
      alpha1g <- alpha1[[g]] # alpha1
    }else{
      alpha1g <- alpha1
    }
    
    pRE <- pREvec[g]
    resid <- Yi - Xi %*% betai - Zi %*% bi
    
    total[(index + 1):(index + pRE)] <- - t(Zi) %*% resid / sigmai
    
    # double check if it is squared
    
    # sum alpha'b
    if (latAsso == "sre") {
      
      sum.alpha1i <- sum.alpha1i + t(alpha1g) %*% bi
      
      dEta1[(index + 1):(index + pRE)] <- as.numeric(alpha1g)
      
    } else if (latAsso == "presentlp") {
      
      Zs_ig <- as.matrix(Zs_i[[g]])
      
      latent <- as.numeric(Zs_ig %*% bi)
      
      sum.alpha1i <- sum.alpha1i + alpha1g * latent
      
      dEta1[(index + 1):(index + pRE)] <- as.numeric(alpha1g * t(Zs_ig))
      
    } else if (latAsso == "present") {
      
      Xs_ig <- as.matrix(Xs_i[[g]])
      Zs_ig <- as.matrix(Zs_i[[g]])
      
      latent <- as.numeric(Xs_ig %*% betai + Zs_ig %*% bi)
      
      sum.alpha1i <- sum.alpha1i + alpha1g * latent
      
      dEta1[(index + 1):(index + pRE)] <- as.numeric(alpha1g * t(Zs_ig))
      
    } else {
      stop("Unknown latent association")
    }
    
    index <- index + pRE
  }
  
  # latent structure for each loop
  latent1 <- as.matrix(sum.alpha1i, nrow = 1)
  CH01 <- as.matrix(CH01)
  
  latent1 <- as.numeric(sum.alpha1i)
  
  eta1 <- as.numeric(Wx %*% gamma1 + latent1)
  
  total <- total +
    as.numeric(CH01 * exp(eta1)) * dEta1 +
    as.numeric(solve(Sig) %*% bfull)
  
  if (Wcmprsk == 1) {
    total <- total - dEta1# adjusts for status == 1
  }
  
  total <- unname(total)
  
  return(total)
  
}