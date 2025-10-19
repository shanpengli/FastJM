Pkmv.us_SF <- function(CIF, data, bl, numBio, pREvec) {
  
  sumalpha1b <- 0
  index <- 0
  
  for (g in 1:numBio) {
    pRE <- pREvec[g]
    if(numBio == 1){ 
      alpha1g <- data$alpha1 
    } else{ 
      alpha1g <- data$alpha1[[g]] 
    }
    
    b_ig <- bl[(index+1):(index+pRE)]
    sumalpha1b <- sumalpha1b + as.numeric(alpha1g%*%b_ig)
    index = index + pRE
  }
  
  CIF/exp(-data$CH01*exp(as.numeric(data$W%*%data$gamma1) + sumalpha1b))

}