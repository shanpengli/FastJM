logLik <- function(data, b) {
  sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*data$sigma)) + data$CH0*exp(data$X2%*%data$gamma + data$nu%*%b) +
    t(b)%*%solve(data$Sig)%*%b
}
