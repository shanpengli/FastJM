logLikCR <- function(data, b) {
  sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*data$sigma)) + 
    data$CH01*exp(data$X2%*%data$gamma1 + data$nu1%*%b) +
    data$CH02*exp(data$X2%*%data$gamma2 + data$nu2%*%b) +
    t(b)%*%solve(data$Sig)%*%b/2
}