Pk.us <- function(CIF, data, bl) {
  CIF/exp(-data$CH01*exp(data$X2%*%data$gamma1 + data$nu1%*%bl)-
            data$CH02*exp(data$X2%*%data$gamma2 + data$nu2%*%bl)
  )
  
}

Pk.us.JMH <- function(CIF, data, bl) {
  p1a <- nrow(data$Sig) - 1
  b <- as.vector(bl[1:p1a])
  w <- bl[p1a+1]
  CIF/exp(-data$CH01*exp(data$X2%*%data$gamma1 + data$alpha1%*%b + data$nu1*w)-
            data$CH02*exp(data$X2%*%data$gamma2 + data$alpha2%*%b + data$nu2*w)
  )
  
}

Pkmv.us <- function(CIF, data, bl, numBio, pREvec) {
  
  sumalpha1b <- 0
  sumalpha2b <- 0
  index <- 0
  
  latAsso <- data$latAsso
  
  for (g in 1:numBio) {
    
    pRE <- pREvec[g]
    
    # random effects for biomarker g
    b_ig <- bl[(index + 1):(index + pRE)]
    
    
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
  
  
  # if(latAsso == "sre"){
  
  # divide CIF by S(s|b)
  CIF/exp(-data$CH01*exp(data$W%*%data$gamma1 + sumalpha1b)-
            data$CH02*exp(data$W%*%data$gamma2 + sumalpha2b)
  )
  # } else {
  #   # CH 0 so denom 1
  #   CIF
  # }
  
  
}
