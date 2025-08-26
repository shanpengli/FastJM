CIF1mv.CR <- function(data, H01, H02, stime, u, bl, numBio, pREvec) {
  a <- nrow(H01)
  b <- nrow(H02)
  CH01 <- cumsum(H01[, 3])
  CH02 <- rep(0, a)
  count <- 1
  i = 1
  while (count <= b & i <= a) {
    if (H02[count, 1] < H01[i, 1]) {
      CH02[i] <- CH02[i] + H02[count, 3]
      count <- count + 1
    } else {
      i <- i + 1
    }
  }
  CH02 <- cumsum(CH02)
  
  CIF1 <- 0
  
  sumalpha1b <- sumalpha2b <- index <- 0
  for(g in 1:numBio){
    pRE <- pREvec[g]
    if(numBio == 1){
      alpha1g <- data$alpha1
      alpha2g <- data$alpha2
    }else{
      alpha1g <- data$alpha1[[g]]
      alpha2g <- data$alpha2[[g]]
    }
    b_ig <- bl[(index+1):(index+pRE)]
    sumalpha1b <- sumalpha1b + alpha1g%*%b_ig # need to double check structure of bl, might be [[i]][index:index+pRE]
    sumalpha2b <- sumalpha2b + alpha2g%*%b_ig
    index = index + pRE
  }
  
  for (i in 1:a) {
     if (stime < H01[i, 1] & u >= H01[i, 1]) {
      if (i >= 2) {
        # CIF1 <- CIF1 + H01[i, 3]*exp(data$W%*%data$gamma1 + data$alpha1%*%bl)*
        #   exp(-CH01[i-1]*exp(data$W%*%data$gamma1 + data$alpha1%*%bl)-
        #         CH02[i-1]*exp(data$W%*%data$gamma2 + data$alpha2%*%bl))
        
        ## will need to change this for fleixble
        
      
        CIF1 <- CIF1 + H01[i, 3]*exp(data$W%*%data$gamma1 + sumalpha1b)*
          exp(-CH01[i-1]*exp(data$W%*%data$gamma1 + sumalpha1b)-
                CH02[i-1]*exp(data$W%*%data$gamma2 + sumalpha2b))
      } else {
        CIF1 <- CIF1 + H01[i, 3]*exp(data$W%*%data$gamma1 + sumalpha1b)
      }
      
    }
  }
  CIF1
}