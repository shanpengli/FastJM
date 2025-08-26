Pkmv.us <- function(CIF, data, bl, numBio, pREvec) {
  
  sumalpha1b <- sumalpha2b <- index <- 0
  
  alpha1 <- alpha2 <- c()
  index = 0
  for(g in 1:numBio){
    pRE <- pREvec[g]
    index = index + pRE
  }
  
  
  index = 0
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
  
  
  CIF/exp(-data$CH01*exp(data$W%*%data$gamma1 + sumalpha1b)-
             data$CH02*exp(data$W%*%data$gamma2 + sumalpha2b)
           )
       
}
