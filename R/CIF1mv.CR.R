CIF1mv.CR <- function(data, H01, H02, stime, u, bl, numBio, pREvec) {
  
  # calculate cum baseline hazard for cause 1
  a <- nrow(H01)
  b <- nrow(H02)
  CH01 <- cumsum(H01[, 3])
  CH02 <- rep(0, a)
  
  # linear scan: sum cause 2 hazard increments between cause 1 event times
  count <- 1
  i = 1
  
  while (count <= b & i <= a) {
    if (H02[count, 1] <= H01[i, 1]) {
      CH02[i] <- CH02[i] + H02[count, 3]
      count <- count + 1
    } else {
      i <- i + 1
    }
  }
  
  CH02 <- cumsum(CH02)
  
  CIF1 <- 0
  sumalpha1b <- 0
  sumalpha2b <- 0
  index <- 0
  
  latAsso <- data$latAsso
  
  # calc association contribution across biomarkers
  for (g in 1:numBio) {
    
    pRE <- pREvec[g]
    
    # random effects for biomarker g
    b_ig <- bl[(index + 1):(index + pRE)]
    
    if (latAsso == "sre") {
      
      alpha1g <- data$alpha1[[g]]
      alpha2g <- data$alpha2[[g]]
      
      latent <- b_ig
      
      sumalpha1b <- sumalpha1b + as.numeric(t(alpha1g) %*% latent)
      sumalpha2b <- sumalpha2b + as.numeric(t(alpha2g) %*% latent)
      
    } else if (latAsso == "presentlp") {
      
      alpha1g <- data$alpha1[g]
      alpha2g <- data$alpha2[g]
      
      Zs_ig <- as.matrix(data$Zs_i[[g]])
      
      # Z_ig(stime) * b_ig
      latent <- as.numeric(Zs_ig %*% b_ig)
      
      sumalpha1b <- sumalpha1b + as.numeric(alpha1g) * latent
      sumalpha2b <- sumalpha2b + as.numeric(alpha2g) * latent
      
    } else if (latAsso == "present") {
      
      alpha1g <- data$alpha1[g]
      alpha2g <- data$alpha2[g]
      
      Xs_ig <- as.matrix(data$Xs_i[[g]])
      Zs_ig <- as.matrix(data$Zs_i[[g]])
      
      
      # X_ig(stime) * beta_g + Z_ig(stime) * b_ig
      latent <- as.numeric(Xs_ig %*% data$beta[[g]] + Zs_ig %*% b_ig)
      
      sumalpha1b <- sumalpha1b + as.numeric(alpha1g) * latent
      sumalpha2b <- sumalpha2b + as.numeric(alpha2g) * latent
      
    } else {
      stop("Unsupported latent association structure: ", latAsso)
    }
    
    # update index for next biomarker
    index <- index + pRE
  }
  
  
  # calculate hazard multipliers
  hazmult1 <- as.numeric(exp(data$W %*% data$gamma1 + sumalpha1b))
  hazmult2 <- as.numeric(exp(data$W %*% data$gamma2 + sumalpha2b))
  
  # approx CIF for cause 1 from pred time to u
  for (i in seq_len(a)) {
    
    if (stime < H01[i, 1] & u >= H01[i, 1]) { # only include events after stime and up to u
      CIF1 <- CIF1 + H01[i, 3] *
        hazmult1 *
        exp(-CH01[i]*hazmult1 -
              CH02[i]*hazmult2)
    }
  }
  
  CIF1
}