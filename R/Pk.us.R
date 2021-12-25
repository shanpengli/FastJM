Pk.us <- function(CIF, data, bl) {
  CIF/exp(-data$CH01*exp(data$X2%*%data$gamma1 + data$nu1%*%bl)-
             data$CH02*exp(data$X2%*%data$gamma2 + data$nu2%*%bl)
           )
       
}