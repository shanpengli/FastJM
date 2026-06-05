getbSig <- function(bSig, data){
  
  
  # sigma <- exp(data$W %*% data$tau + w)  # calculates sigma; don't need assume homogeneous error variance
  # sum over G
  # Sig: covariance matrix of B
  # sigma: vector of error variance for each biomarker
  # unique to each value

  latAsso <- data$latAsso
  s <- data$s
  
  Y <- data$Y
  X <- data$X # update so both biomarkers accountted for
  Z <- data$Z
  
  Xs_i <- data$Xs_i
  Zs_i <- data$Zs_i
  
  beta <- data$beta
  
  alphaList <- data$alphaList
  sigma <- data$sigma
  Sig <- data$Sig
  
  # Sig11 <- SigList[[1]]
  # Sig22<- SigList[[2]] # first biomarker/ #second biomarker
  
  # survival
  CH01 <- data$CH01
  CH02 <- data$CH02
  

  HAZ01 <- data$HAZ01
  HAZ02 <- data$HAZ02
  Wcmprsk <- data$Wcmprsk
  Wx <- as.matrix(data$W)
  if(ncol(Wx) ==1){
    Wx <- t(Wx)
  }
  gamma1 <- as.matrix(data$gamma1) # vector
  gamma2 <- as.matrix(data$gamma2)
  
  # turn into list/loop
  pREvec <- c()
  if(is.list(Z)){
    for(g in 1:length(Z)){
      pREvec[g] <- ncol(Z[[g]])
    }
  }else{
    g <- 1
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
  sum.alpha2i <- 0
  
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
    alpha1 <- alphaList[[1]] # risk 1
    alpha2 <- alphaList[[2]] # risk 2
    # gets for each biomarker
    
    if(is.list(alpha1)){
      alpha1g <- alpha1[[g]] # alpha1
      alpha2g <- alpha2[[g]] # alpha2
    }else{
      alpha1g <- alpha1
      alpha2g <- alpha2
    }
    
    # log likelihood part 1
    total <- total + sum((Yi - Xi %*% betai - Zi %*% bi)^2 / (2 * sigmai) + 0.5 * log(sigmai))
    
    
    # sum alpha'b
    
    if(latAsso == "sre"){
      sum.alpha1i <- sum.alpha1i + t(alpha1g) %*% bi #alpha1
      sum.alpha2i <- sum.alpha2i + t(alpha2g) %*% bi #alpha2
    } else if(latAsso  == "presentlp"){
        Zs_ig <- as.matrix(Zs_i[[g]])
          # if(length(alpha1g)== 1){
        latent <- as.numeric(Zs_ig %*% bi)
            sum.alpha1i <- sum.alpha1i + alpha1g * latent #alpha1
            sum.alpha2i <- sum.alpha2i + alpha2g * latent #alpha2
    } else if(latAsso == "present"){
      Xs_ig <- as.matrix(Xs_i[[g]])
      Zs_ig <- as.matrix(Zs_i[[g]])
      
      latent <- as.numeric(Xs_ig %*% betai + Zs_ig %*% bi)
        sum.alpha1i <- sum.alpha1i + alpha1g * latent #alpha1
        sum.alpha2i <- sum.alpha2i + alpha2g * latent #alpha2
    } 
  }


  # latent structure for each loop
  latent1 <- as.numeric(sum.alpha1i)
  latent2 <- as.numeric(sum.alpha2i)
  
  total <- total + CH01 * exp(Wx %*% gamma1 + latent1) + ## part 2 change this part
    CH02 * exp(Wx %*% gamma2 + latent2) +
    0.5 * q *log(2*pi) +
    0.5*log(det(Sig)) + t(bfull) %*% solve(Sig) %*% bfull / 2  # part 3
  
  # print(total)
  
  if (Wcmprsk == 1) {
    total <- total - log(HAZ01) - (Wx %*% gamma1 + latent1) # adjusts for status == 1
  }
  
  if (Wcmprsk == 2) {
    total <- total - log(HAZ02) - (Wx %*% gamma2 + latent2)  # adjusts for status == 2
  }

  
  total <- unname(total)
  
  return(total)
  
}