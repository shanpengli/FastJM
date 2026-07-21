logLik <- function(data, b) {
  sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*data$sigma)) + data$CH0*exp(data$X2%*%data$gamma + data$nu%*%b) +
    t(b)%*%solve(data$Sig)%*%b/2
}

logLik.learn <- function(data, b) {
  total <- sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*data$sigma) + 0.5*log(data$sigma)) + 
    data$CH01*exp(data$X2%*%data$gamma1 + data$nu1%*%b) +
    t(b)%*%solve(data$Sig)%*%b/2 + 0.5*log(det(data$Sig))
  if (data$D == 1) {
    total <- total - log(data$HAZ01) - (data$X2%*%data$gamma1 + data$nu1%*%b)
    total
  } else {
    total
  }
}

logLik.JMH <- function(data, bw) {
  p1a <- nrow(data$Sig) - 1
  b <- as.vector(bw[1:p1a])
  w <- bw[p1a+1]
  sigma <- exp(data$W%*%data$tau + w)
  sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*sigma) + 0.5*log(sigma)) + data$CH0*exp(data$X2%*%data$gamma + data$alpha%*%b + data$nu*w) +
    t(bw)%*%solve(data$Sig)%*%bw/2 + 0.5*log(det(data$Sig))
}

logLik.learn.JMH <- function(data, bw) {
  p1a <- nrow(data$Sig) - 1
  b <- as.vector(bw[1:p1a])
  w <- bw[p1a+1]
  sigma <- exp(data$W%*%data$tau + w)
  total <- sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*sigma) + 0.5*log(sigma)) + 
    data$CH01*exp(data$X2%*%data$gamma1 + data$alpha1%*%b + data$nu1*w) +
    t(bw)%*%solve(data$Sig)%*%bw/2 + 0.5*log(det(data$Sig))
  if (data$D == 1) {
    total <- total - log(data$HAZ01) - (data$X2%*%data$gamma1 + data$alpha1%*%b + data$nu1*w)
    total
  } else {
    total
  }
}

logLik.learn.mv <- function(data, bSig) {
  
  Y <- data$Y
  X <- data$X # update so both biomarkers accountted for
  Z <- data$Z
  
  alphaList <- data$alphaList
  sigma <- data$sigma
  Sig <- data$Sig
  beta <- data$beta
  
  # Sig11 <- SigList[[1]]
  # Sig22<- SigList[[2]] # first biomarker/ #second biomarker
  
  # survival
  CH01 <- data$CH01
  
  HAZ01 <- data$HAZ01
  Wcmprsk <- data$Wcmprsk
  Wx <- data$W
  if(ncol(Wx) ==1){
    Wx <- t(Wx)
  }
  gamma1 <- as.matrix(data$gamma1) # vector
  
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
    # gets for each biomarker
    
    if(is.list(alpha1)){
      alpha1g <- alpha1[[g]] # alpha1
    }else{
      alpha1g <- alpha1
    }
    
    # log likelihood part 1
    total <- total + sum((Yi - Xi %*% betai - Zi %*% bi)^2 / (2 * sigmai) + 0.5 * log(sigmai))
    
    
    # sum alpha'b
    sum.alpha1i <- sum.alpha1i + t(alpha1g) %*% bi #alpha1
  }
  
  # latent structure for each loop
  latent1 <- as.numeric(sum.alpha1i)
  
  total <- total + CH01 * exp(Wx %*% gamma1 + latent1) + ## part 2 change this part
    0.5 * q *log(2*pi) +
    0.5*log(det(Sig)) + t(bfull) %*% solve(Sig) %*% bfull / 2  # part 3
  # print(total)
  if (Wcmprsk == 1) {
    total <- total - log(HAZ01) - (Wx %*% gamma1 + latent1) # adjusts for status == 1
  }
  total <- unname(total)
  
  return(total)
}
