getbSigSF <- function(bSig, data){

  # sigma <- exp(data$W %*% data$tau + w)  # calculates sigma; don't need assume homogeneous error variance
  # sum over G
  # Sig: covariance matrix of B
  # sigma: vector of error variance for each biomarker
  # unique to each value
  
  names(data) <- c("beta", "gamma1", "alphaList",
                   "sigma", "Z", "X", "Y", "Sig", # "b", "Sig",
                   "CH01",
                   "HAZ01", "mdata", "mdataS", "Wcmprsk", "Wx")
  
  # don't need mdata
  Y <- data$Y
  X <- data$X # update so both biomarkers accountted for
  Z <- data$Z
  
  beta <- data$beta
  
  alphaList <- data$alphaList
  sigma <- data$sigma
  # SigList <- data$SigList # need to be 4x4
  
  Sig <- data$Sig
  
  # survival
  CH01 <- data$CH01
  HAZ01 <- data$HAZ01
  mdata <- data$mdata
  mdataS <- data$mdataS
  Wcmprsk <- data$Wcmprsk
  Wx <- as.matrix(data$Wx)
  gamma1 <- as.matrix(data$gamma1) # vector
  
  # turn into list/loop
  pREvec <- c()
  if(is.list(Z)){
    for(g in 1:length(Z)){
      pREvec[g] <- ncol(Z[[g]])
    }
  }else{
    pREvec[g] <- ncol(Z)
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
  
  total <- 0
  sum.alpha1i <- 0
  
  # need to generalize here
  
  # longitudinal portion
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
    mdatag <- mdata[[g]]
    mdataSg <- mdataS[[g]]
    alpha1 <- alphaList[[1]] # risk 1
    # gets for each biomarker
    
    if(is.list(alpha1)){
      alpha1g <- alpha1[[g]] # alpha1
    }else{
      alpha1g <- alpha1
    }
    
    # log likelihood part 1
    total <- total + sum((Yi - Xi %*% betai - Zi %*% bi)^2 / (2 * sigmai) + 0.5 * log(sigmai))
    
    # double check if it is squared
    
    # sum alpha'b
    sum.alpha1i <- sum.alpha1i + t(alpha1g) %*% bi #alpha1
  }
  
  # latent structure for each loop
  latent1 <- as.matrix(sum.alpha1i, nrow = 1)
  CH01 <- as.matrix(CH01)
  
  # CH01 Might be wrong here
  
  total <- total + CH01 * exp(Wx%*% gamma1 + latent1) + ## part 2 change this part
    0.5 * q *log(2*pi) +
    0.5*log(det(Sig)) + t(bfull) %*% solve(Sig) %*% bfull / 2  # part 3
  
  if (Wcmprsk == 1) {
    total <- total - log(HAZ01) - (Wx %*% gamma1 + latent1) # adjusts for status == 1
  }
  
  total <- unname(total)
  
  return(total)
  
}