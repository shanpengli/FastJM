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
  
  
  sumalpha1b <- sumalpha2b <- index <- 0
  for(g in 1:numBio){
    pRE <- pREvec[g]
    if(numBio == 1){
      alpha1g <- data$alpha[[1]]
      alpha2g <- data$alpha[[2]]
    }else{
      alpha1g <- data$alpha[[1]][[g]]
      alpha2g <- data$alpha[[2]][[g]]
    }
    b_ig <- bl[(index+1):(index+pRE)]
    sumalpha1b <- sumalpha1b + alpha1g%*%b_ig # need to double check structure of bl, might be [[i]][index:index+pRE]
    sumalpha2b <- sumalpha2b + alpha2g%*%b_ig
    index = index + pRE
  }
  
  for (i in 1:b) {
    if (stime <= H02[i, 1] & u >= H02[i, 1]) {
      if (i >= 2) {
        
        CIF2 <- CIF2 + H02[i, 3]*exp(data$W%*%data$gamma2 + sumalpha2b)*
          exp(-CH01[i-1]*exp(data$W%*%data$gamma1 + sumalpha1b)-
                CH02[i-1]*exp(data$W%*%data$gamma2 + sumalpha2b))
        
        
      } else {
        CIF2 <- CIF2 + H02[i, 3]*exp(data$W%*%data$gamma2 + sumalpha2b)
      }
      
      
    }
  }
  CIF2
}