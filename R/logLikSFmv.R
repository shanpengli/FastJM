logLikSFmv <- function(data, b) {
  
  total <- 0
  index <- 1
  numBio <- length(data$YList)
  sum.alpha1i <- 0
  
  for(g in 1:numBio){
    
    bg <- b[index:(index + length(data$alpha1[[g]])-1)]
    
    total = total + 
      sum(((data$YList[[g]]-data$XList[[g]]%*%data$beta[[g]]) - data$ZList[[g]]%*%bg)^2/(2*as.numeric(data$sigma[g])))
    index = index + length(data$alpha1[[g]])
    
    if(is.list(data$alpha1)){
      alpha1g <- data$alpha1[[g]] # alpha1 bio g
    }else{
      alpha1g <- data$alpha1
    }
    
    sum.alpha1i <- sum.alpha1i + t(alpha1g) %*% bg #alpha1
  }
  
  total = total + data$CH01*exp(data$W%*%data$gamma1 + sum.alpha1i)

 
  total = total  + t(b)%*%solve(data$Sig)%*%b/2
  
  return(unname(total))

}