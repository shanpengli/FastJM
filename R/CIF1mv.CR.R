CIF1mv.CR <- function(data, H01, H02, s, u, bl, numBio, pREvec) {
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
  
  for (i in 1:a) {
    if (s < H01[i, 1] & u >= H01[i, 1]) {
      if (i >= 2) {
        # CIF1 <- CIF1 + H01[i, 3]*exp(data$W%*%data$gamma1 + data$alpha1%*%bl)*
        #   exp(-CH01[i-1]*exp(data$W%*%data$gamma1 + data$alpha1%*%bl)-
        #         CH02[i-1]*exp(data$W%*%data$gamma2 + data$alpha2%*%bl))
        
        ## will need to change this for fleixble
        
        
        sumalpha1b <- sumalpha2b <- index <- 0
        for(g in 1:numBio){
          pRE <- pREvec(g)
          b_ig <- bl[[i]][[g]]
          sumalpha1b <- sumalpha1b + alpha1[index:index+pRE]%*%b_ig # need to double check structure of bl, might be [[i]][index:index+pRE]
          sumalpha2b <- sumalpha2b + alpha2[index:index+pRE]%*%b_ig
          index = index + pRE
        }
        
        CIF1 <- CIF1 + H01[i, 3]*exp(data$W%*%data$gamma1 + sumalpha1b)*
          exp(-CH01[i-1]*exp(data$W%*%data$gamma1 + sumalpha1b)-
                CH02[i-1]*exp(data$W%*%data$gamma2 + sumalpha2b))
      } else {
        CIF1 <- CIF1 + H01[i, 3]*exp(data$W%*%data$gamma1 + sumalpha1b)
      }
      
    } else next
  }
  CIF1
}