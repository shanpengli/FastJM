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