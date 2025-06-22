CIF2.CR <- function(data, H01, H02, s, u, bl) {
  a <- nrow(H01)
  b <- nrow(H02)
  CH02 <- cumsum(H02[, 3])
  CH01 <- rep(0, b)
  count <- 1
  i = 1
  while (count <= a & i <= b) {
    if (H01[count, 1] <= H02[i, 1]) {
      CH01[i] <- CH01[i] + H01[count, 3]
      count <- count + 1
    } else {
      i <- i + 1
    }
  }
  CH01 <- cumsum(CH01)
  
  CIF2 <- 0
  for (i in 1:b) {
    if (s <= H02[i, 1] & u >= H02[i, 1]) {
      if (i >= 2) {
        CIF2 <- CIF2 + H02[i, 3]*exp(data$X2%*%data$gamma2 + data$nu2%*%bl)*
          exp(-CH01[i-1]*exp(data$X2%*%data$gamma1 + data$nu1%*%bl)-
                CH02[i-1]*exp(data$X2%*%data$gamma2 + data$nu2%*%bl))
      } else {
        CIF2 <- CIF2 + H02[i, 3]*exp(data$X2%*%data$gamma2 + data$nu2%*%bl)
      }

    }
  }
  CIF2
}