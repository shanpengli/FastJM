CIF1mv.SF <- function(data, H01, stime, u, bl, numBio, pREvec) {
  
  a <- nrow(H01)
  CH01 <- cumsum(H01[, 3])  # cumulative baseline hazard for event 1
  
  CIF1 <- 0
  sumalpha1b <- 0
  index <- 0
  
  for (g in 1:numBio) {
    pRE <- pREvec[g]
    if(numBio == 1){
      alpha1g <- data$alpha1
    }else{
      alpha1g <- data$alpha1[[g]]
    }
    b_ig <- bl[(index + 1):(index + pRE)]
    sumalpha1b <- sumalpha1b + as.numeric(alpha1g %*% b_ig)
    index <- index + pRE
  }
  
  wgamma <- as.numeric(data$W%*%data$gamma1) 
  for (i in 1:a) { 
    if (stime < H01[i, 1] && u >= H01[i, 1]) { 
      if (i >= 2) { 
        CIF1 <- CIF1 + H01[i, 3]*exp(wgamma + sumalpha1b)* exp(-CH01[i-1]*exp(wgamma + sumalpha1b)) 
      } else { 
          CIF1 <- CIF1 + H01[i, 3]*exp(wgamma + sumalpha1b) 
      } 
    } 
  }
  
  CIF1
}



