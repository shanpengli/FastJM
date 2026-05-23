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
  
  # Splitting the stacked random-effects vector into 
  # biomarker-specific pieces
  b_list <- split(as.vector(bl), rep(seq_len(numBio), times = pREvec))
  
  # For saving sum_g alpha_g^T b_g
  alpha_b_terms <- numeric(numBio)
  
  # Looping through biomarkers
  for (g in seq_len(numBio)) {
    alpha_b_terms[g] <- data$alpha1[[g]] %*% b_list[[g]]
  }
  
  alpha_b_sum <- sum(alpha_b_terms)
  
  (exp(-data$CH01 * exp(data$W %*% data$gamma1 + alpha_b_sum)) -
      exp(-CH0u * exp(data$W %*% data$gamma1 + alpha_b_sum))) /
    exp(-data$CH01 * exp(data$W %*% data$gamma1 + alpha_b_sum))
}