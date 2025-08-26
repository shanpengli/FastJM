logLikCRmv <- function(data, b) {
  
  total <- 0
  index <- 1
  numBio <- length(data$Y)
  
  
  sum.alpha1i <- 0
  sum.alpha2i <- 0
  

  
  for(g in 1:numBio){
    
    bg <- b[index:(index + length(data$alpha1[[g]])-1)]
    
    total = total + 
      sum(((data$YList[[g]]-data$XList[[g]]%*%data$beta[[g]]) - data$ZList[[g]]%*%bg)^2/(2*data$sigma[g]))
    index = index + length(data$alpha1[[g]])
    
    if(is.list(data$alpha1)){
      alpha1g <- data$alpha1[[g]] # alpha1 bio g
      alpha2g <- data$alpha2[[g]] # alpha2 bio g
    }else{
      alpha1g <- data$alpha1
      alpha2g <- data$alpha2
    }
    
    sum.alpha1i <- sum.alpha1i + t(alpha1g) %*% bg #alpha1
    sum.alpha2i <- sum.alpha2i + t(alpha2g) %*% bg #alpha2
  }
  
  total = total + data$CH01*exp(data$W%*%data$gamma1 + sum.alpha1i) +
    data$CH02*exp(data$W%*%data$gamma2 + sum.alpha2i)

 
  total = total  + t(b)%*%solve(data$Sig)%*%b/2
  
  return(unname(total))

}