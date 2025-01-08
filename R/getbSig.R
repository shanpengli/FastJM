getbSig <- function(bSig, data){
  
  
  # sigma <- exp(data$W %*% data$tau + w)  # calculates sigma; don't need assume homogeneous error variance
  # sum over G
  # Sig: covariance matrix of B
  # sigma: vector of error variance for each biomarker
  # unique to each value
  
  names(data) <- c("beta", "gamma1", "gamma2", "alphaList",
                   "sigma", "Z", "X", "Y", "Sig", # "b", "Sig",
                   "CH01", "CH02",
                   "HAZ01", "HAZ02", "mdata", "mdataS", "Wcmprsk", "Wx")
  
  
  
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
  CH02 <- data$CH02
  HAZ01 <- data$HAZ01
  HAZ02 <- data$HAZ02
  mdata <- data$mdata
  mdataS <- data$mdataS
  Wcmprsk <- data$Wcmprsk
  Wx <- t(as.matrix(data$Wx))
  gamma1 <- as.matrix(data$gamma1) # vector
  gamma2 <- as.matrix(data$gamma2)
  
  # p12 <- nrow(Z[[1]])
  # p22<- nrow(Z[[2]])
  
  
  # turn into list/loop
  p1a <- ncol(Z[[1]]) # dim of random effect
  p2a <- ncol(Z[[2]])
  
  
  # for(i in length(Z))P{}
  # numSub <- length(mdata[[1]])
  
  q = 0
  
  for(i in 1:length(Z)){
    q = q + ncol(Z[[i]])
  }
  
  
  # b <- list(
  #   b1 = matrix(bSig[1:(2*p12)], nrow = p12, ncol = p1a),
  #   b2 = matrix(bSig[(2*p12+1):(2*p12+2*p22)], nrow = p22, ncol = p1a)
  # )
  #
  # Sig <- matrix(bSig[(2*p12+2*p22+1):length(bSig)], nrow = p1a, ncol = p1a
  
  b <- list(
    b1 = matrix(bSig[1:p1a], nrow = p1a, ncol = 1),
    b2 = matrix(bSig[(p1a+1):(p1a+p2a)], nrow = p2a, ncol = 1)
  )
  
  
  
  #Sig <- matrix(bSig[(2*q+2*q+1):length(bSig)], nrow = p1a, ncol = p1a)
  
  total <- 0
  sum.alpha1i <- 0
  sum.alpha2i <- 0
  
  # need to generalize here
  bfull <- matrix(ncol = 1)
  
  
  # longitudinal portion
  for (g in 1:length(Y)) {
    
    Yi <- as.matrix(Y[[g]])
    Xi <- as.matrix(X[[g]])
    betai <- as.matrix(beta[[g]])
    Zi <- as.matrix(Z[[g]])
    bi <- as.matrix(b[[g]])
    sigmai <- sigma[[g]]
    mdatag <- mdata[[g]]
    mdataSg <- mdataS[[g]]
    alpha1 <- alphaList[[1]] # risk 1
    alpha2 <- alphaList[[2]] # risk 2
    # gets for each biomarker
    alpha1g <- alpha1[[g]] # bio1
    alpha2g <- alpha2[[g]] # bio2
    
    # calculate zbig
    # need to adjust for repeated measures
    # Zbig <- matrix(NA, nrow = nrow(Zi), ncol = 1)
    # for (i in 1:length(mdatag)) {
    #   repM <- mdatag[i]
    #   for (rep in 1:repM) {
    #     # each row * each bi
    #     Zbig[mdataSg[i] + rep - 1, ] <- t(Zi[mdataSg[i] + rep - 1, ]) %*% bi[i, ]
    #   }
    # }
    
    # log likelihood part 1
    total <- total + sum((Yi - Xi %*% betai - Zi %*% bi)^2 / (2 * sigmai) + 0.5 * log(sigmai))
    
    # double check if it is squared
    
    # sum alpha'b
    sum.alpha1i <- sum.alpha1i + t(alpha1g) %*% bi #bio 1
    sum.alpha2i <- sum.alpha2i + t(alpha2g) %*% bi #bio 2
    
    # get the full re vector
    bfull <- rbind(bfull, bi)
  }
  
  # remove NA row
  bfull<- t(bfull[-1,])
  
  # latent structure for each loop
  latent1 <- as.matrix(rep(sum.alpha1i, nrow(Wx)), nrow = 1, ncol = nrow(Wx))
  latent2 <- as.matrix(rep(sum.alpha2i, nrow(Wx)), nrow = 1, ncol = nrow(Wx))
  CH01 <- as.matrix(CH01)
  CH02 <- as.matrix(CH02)
  
  
  
  # CH01 Might be wrong here
  
  
  total <- total +CH01 * exp(Wx%*% gamma1 + latent1) + ## part 2 change this part
    CH02 * exp(Wx %*% gamma2 + latent2) +
    0.5 * q *log(2*pi) +
    0.5*log(det(Sig)) + bfull %*% solve(Sig) %*% t(bfull) / 2  # part 3
  # need to update SIG/fix dimension
  # need to account for
  
  if (Wcmprsk == 1) {
    # should be HAZ?, double check though
    # total <- total - log(HAZ01[i]) - (W %*% gamma1 + sum.alpha1i)[1,1]  # adjusts for status == 1
    total <- total - log(HAZ01) - (Wx %*% gamma1 + latent1) # adjusts for status == 1
  }
  
  if (Wcmprsk == 2) {
    # total <- total - log(HAZ02[i]) - (W %*% gamma2 + sum.alpha2i)[1,1]  # adjusts for status == 2
    total <- total - log(HAZ02) - (Wx %*% gamma2 + latent2)  # adjusts for status == 2
  }
  
  
  
  
  total <- unname(total)
  return(total)
  
  # for(k in 1:2){
  #   if (data$status == k) {
  #     total <- total - log(data$HAZ0[k]) - (W %*% gammag + alphag %*% b)  # adjusts for status == 1
  #   }
  # }
  
}