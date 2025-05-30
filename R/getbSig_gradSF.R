getbSig_gradSF <- function(bSig, data){
  
  
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
  
  # how to get Sig1,2
  # can do (0,0,0,0) for now; email shangpeng later
  
  Sig <- data$Sig
  
  # Sig11 <- SigList[[1]]
  # Sig22<- SigList[[2]] # first biomarker/ #second biomarker
  
  # survival
  CH01 <- data$CH01
  HAZ01 <- data$HAZ01
  mdata <- data$mdata
  mdataS <- data$mdataS
  Wcmprsk <- data$Wcmprsk
  Wx <- as.matrix(data$Wx)
  gamma1 <- as.matrix(data$gamma1) # vector
  
  # p12 <- nrow(Z[[1]])
  # p22<- nrow(Z[[2]])
  
  
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
  
  total <- c()
  sum.alpha1i <- 0
  sum.alpha2i <- 0
  
  # need to generalize here
  
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
    mdatag <- mdata[[g]]
    mdataSg <- mdataS[[g]]
    alpha1 <- alphaList[[1]] # risk 1
    # gets for each biomarker
    
    if(is.list(alpha1)){
      alpha1g <- alpha1[[g]] # alpha1
    }else{
      alpha1g <- alpha1
    }
    
    pRE <- pREvec[g]
    total[(index+1):(index+pRE)] <- - 2*t(Zi) %*% (Yi - Xi %*% betai - Zi %*% bi) / (2 * sigmai)
    index <- index + pRE
    
    # double check if it is squared
    
    # sum alpha'b
    sum.alpha1i <- sum.alpha1i + t(alpha1g) %*% bi #alpha1
  }
  
  # latent structure for each loop
  latent1 <- sum.alpha1i
  CH01 <- as.matrix(CH01)
  
  # CH01 Might be wrong here
  
  total <- total + as.numeric(CH01 * exp(Wx%*% gamma1 + latent1))*unlist(alpha1) + ## part 2 change this part
 + solve(Sig) %*% bfull  # part 3
  
  if (Wcmprsk == 1) {
    total <- total - unlist(alpha1) # adjusts for status == 1
  }
  
  total <- unname(total)
  
  return(total)
  
}