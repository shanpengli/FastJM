logLikCR <- function(data, b) {
  sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*data$sigma)) + 
    data$CH01*exp(data$X2%*%data$gamma1 + data$nu1%*%b) +
    data$CH02*exp(data$X2%*%data$gamma2 + data$nu2%*%b) +
    t(b)%*%solve(data$Sig)%*%b/2
}

logLikCR.learn <- function(data, b) {
  total <- sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*data$sigma) + 0.5*log(data$sigma)) + 
    data$CH01*exp(data$X2%*%data$gamma1 + data$nu1%*%b) +
    data$CH02*exp(data$X2%*%data$gamma2 + data$nu2%*%b) +
    t(b)%*%solve(data$Sig)%*%b/2 + 0.5*log(det(data$Sig))
  if (data$D == 1) {
    total <- total - log(data$HAZ01) - (data$X2%*%data$gamma1 + data$nu1%*%b)
    total
  } else if (data$D == 2) {
    total <- total - log(data$HAZ02) - (data$X2%*%data$gamma2 + data$nu2%*%b)
    total
  } else {
    total
  }
}