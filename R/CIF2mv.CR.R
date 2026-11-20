CIF2mv.CR <- function(data, H01, H02, stime, u, bl, numBio, pREvec) {
  a <- nrow(H01)
  b <- nrow(H02)
  CH02 <- cumsum(H02[, 3])
  CH01 <- rep(0, b)
  
  count <- 1
  i = 1
  while (count <= a & i <= b) {
    if (H01[count, 1] <= H02[i, 1]) {
      CH01[i] <- CH01[i] + H01[count, 3]
      count <- count + 1
    } else {
      i <- i + 1
    }
  }
  CH01 <- cumsum(CH01)
  
  CIF2 <- 0
  
  sumalpha1b <- 0
  sumalpha2b <- 0
  index <- 0
  
  latAsso <- data$latAsso
  
  for (g in 1:numBio) {
    
    pRE <- pREvec[g]
    
    # random effects for biomarker g
    b_ig <- bl[(index + 1):(index + pRE)]
    
    # association contribution
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
      
      latent <- as.numeric(Zs_ig %*% b_ig)
      
      sumalpha1b <- sumalpha1b + as.numeric(alpha1g) * latent
      sumalpha2b <- sumalpha2b + as.numeric(alpha2g) * latent
      
    } else if (latAsso == "present") {
      
      alpha1g <- data$alpha1[g]
      alpha2g <- data$alpha2[g]
      
      Xs_ig <- as.matrix(data$Xs_i[[g]])
      Zs_ig <- as.matrix(data$Zs_i[[g]])
      
      latent <- as.numeric(Xs_ig %*% data$beta[[g]] + Zs_ig %*% b_ig)
      
      sumalpha1b <- sumalpha1b + as.numeric(alpha1g) * latent
      sumalpha2b <- sumalpha2b + as.numeric(alpha2g) * latent
      
    } else {
      stop("Unsupported latent association structure: ", latAsso)
    }
    
    index <- index + pRE
  }
  
  
  hazmult1 <- as.numeric(exp(data$W %*% data$gamma1 + sumalpha1b))
  hazmult2 <- as.numeric(exp(data$W %*% data$gamma2 + sumalpha2b))
  
  for (i in seq_len(b)) {
    if (stime < H02[i, 1] & u >= H02[i, 1]) {
      CIF2 <- CIF2 + H02[i, 3] *
        hazmult2 *
        exp(-CH01[i]*hazmult1 -
              CH02[i]*hazmult2)
    }
  }
  
  CIF2
}