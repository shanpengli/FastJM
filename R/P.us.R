P.us <- function(data, CH0u, bl) {
  (exp(-data$CH0*exp(data$X2%*%data$gamma + data$nu%*%bl)) - 
     exp(-CH0u*exp(data$X2%*%data$gamma + data$nu%*%bl)))/
    exp(-data$CH0*exp(data$X2%*%data$gamma + data$nu%*%bl))
}

P.us.JMH <- function(data, CH0u, bwl) {
  p1a <- nrow(data$Sig) - 1
  bl <- as.vector(bwl[1:p1a])
  wl <- bwl[p1a+1]
  (exp(-data$CH0*exp(data$X2%*%data$gamma + data$alpha%*%bl + data$nu*wl)) - 
      exp(-CH0u*exp(data$X2%*%data$gamma + data$alpha%*%bl + data$nu*wl)))/
    exp(-data$CH0*exp(data$X2%*%data$gamma + data$alpha%*%bl + data$nu*wl))
}

P.us.mvjmcs <- function(data, CH0u, bl, numBio, pREvec) {
  
  ## EDIT THIS PART
  
  latAsso <- data$latAsso
  
  # split stacked RE by biomarker
  b_list <- split(as.vector(bl), rep(seq_len(numBio), times = pREvec))
  
  
  # Calculate association term
  alpha_b_sum <- 0
  
  for (g in seq_len(numBio)) {
    
    bg <- b_list[[g]]
    
    if (latAsso == "sre") {
      
      # Shared random effects
      alpha_b_sum <- alpha_b_sum +
        as.numeric(data$alpha1[[g]] %*% bg)
      
    } else if (latAsso == "presentlp") {
      
      # Present latent process at landmark s
      latent <- as.numeric(data$Zs_i[[g]] %*% bg)
      
      alpha_b_sum <- alpha_b_sum +
        data$alpha1[[g]] * latent
      
    } else if (latAsso == "present") {
      
      # Present value at landmark s
      latent <- as.numeric(
        data$Xs_i[[g]] %*% data$beta[[g]] +
          data$Zs_i[[g]] %*% bg
      )
      
      alpha_b_sum <- alpha_b_sum +
        data$alpha1[[g]] * latent
      
    } else {
      stop("Unsupported association structure: ", latAsso)
    }
  }
  
  # if(latAsso == "sre"){
  #   (exp(-data$CH01 * exp(data$W %*% data$gamma1 + alpha_b_sum)) -
  #      exp(-CH0u * exp(data$W %*% data$gamma1 + alpha_b_sum))) /
  #     exp(-data$CH01 * exp(data$W %*% data$gamma1 + alpha_b_sum))
  # }
  # # if use landmarking
  # else{
  #   1-exp(-CH0u * exp(data$W %*% data$gamma1 + alpha_b_sum))
  # }
  
  # calc hazard multiplier
  hazmult1 <- exp(data$W %*% data$gamma1 + alpha_b_sum)
  
  # calc conditional failure probability from prediction time to u
  1 - exp(-(CH0u - data$CH01) * hazmult1)
  
}