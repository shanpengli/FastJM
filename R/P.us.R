P.us <- function(data, CH0u, bl) {
  (exp(-data$CH0*exp(data$X2%*%data$gamma + data$nu%*%bl)) - 
     exp(-CH0u*exp(data$X2%*%data$gamma + data$nu%*%bl)))/
    exp(-data$CH0*exp(data$X2%*%data$gamma + data$nu%*%bl))
}