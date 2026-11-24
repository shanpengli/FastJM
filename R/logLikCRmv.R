logLikCRmv <- function(data, b) {
  
  total <- 0
  index <- 1
  numBio <- length(data$YList)
  
  sum.alpha1i <- 0
  sum.alpha2i <- 0
  
  latAsso <- data$latAsso
  
  # only needed for present/presentlp
  if (latAsso %in% c("present", "presentlp")) {
    Xs_i <- data$Xs_i
    Zs_i <- data$Zs_i
  }
  
  for(g in 1:numBio){
    
    # number of random effects for biomarker g
    pREg <- ncol(data$ZList[[g]])
    
    # random effects for biomarker g
    bg <- b[index:(index + pREg - 1)]
    
    # longitudinal contribution
    total = total + 
      sum(((data$YList[[g]]-data$XList[[g]]%*%data$beta[[g]]) - data$ZList[[g]]%*%bg)^2/(2*data$sigma[g]))
    
    # association contribution
    if (latAsso == "sre") {
      
      alpha1g <- data$alpha1[[g]]
      alpha2g <- data$alpha2[[g]]
      
      latent <- bg
      
    } else if (latAsso == "presentlp") {
      
      alpha1g <- data$alpha1[g]
      alpha2g <- data$alpha2[g]
      
      Zs_ig <- as.matrix(Zs_i[[g]])
      
      latent <- as.numeric(Zs_ig %*% bg)
      
    } else if (latAsso == "present") {
      
      alpha1g <- data$alpha1[g]
      alpha2g <- data$alpha2[g]
      
      Xs_ig <- as.matrix(Xs_i[[g]])
      Zs_ig <- as.matrix(Zs_i[[g]])
      
      latent <- as.numeric(Xs_ig %*% data$beta[[g]] + Zs_ig %*% bg)
      
    }  else {
      stop("Unsupported latent association structure: ", latAsso)
    }
    
    sum.alpha1i <- sum.alpha1i + t(alpha1g) %*% latent # alpha1
    sum.alpha2i <- sum.alpha2i + t(alpha2g) %*% latent # alpha2
    
    index = index + pREg
  }
  
  latent1 <- as.numeric(sum.alpha1i)
  latent2 <- as.numeric(sum.alpha2i)
  
  # association contribution
  if(latAsso == "sre"){
    # only consider for sre; CH01 is cum hazard at time = s--> CH01 = 0 ow
    total <- total + data$CH01*exp(data$W%*%data$gamma1 + latent1) +
      data$CH02*exp(data$W%*%data$gamma2 + latent2)
  }
  
  # RE contribution
  total = total  + t(b)%*%solve(data$Sig)%*%b/2
  
  return(unname(total))
  
}