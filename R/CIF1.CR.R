CIF1.CR <- function(data, H01, H02, s, u, bl) {
  a <- nrow(H01)
  b <- nrow(H02)
  CH01 <- cumsum(H01[, 3])
  CH02 <- rep(0, a)
  count <- 1
  i = 1
  while (count <= b & i <= a) {
    if (H02[count, 1] < H01[i, 1]) {
      CH02[i] <- CH02[i] + H02[count, 3]
      count <- count + 1
    } else {
      i <- i + 1
    }
  }
  CH02 <- cumsum(CH02)
  
  CIF1 <- 0
  for (i in 1:a) {
    if (s < H01[i, 1] & u >= H01[i, 1]) {
      if (i >= 2) {
        CIF1 <- CIF1 + H01[i, 3]*exp(data$X2%*%data$gamma1 + data$nu1%*%bl)*
          exp(-CH01[i-1]*exp(data$X2%*%data$gamma1 + data$nu1%*%bl)-
                CH02[i-1]*exp(data$X2%*%data$gamma2 + data$nu2%*%bl))
      } else {
        CIF1 <- CIF1 + H01[i, 3]*exp(data$X2%*%data$gamma1 + data$nu1%*%bl)
      }
      
    }
  }
  CIF1
}

CIF1.CR.JMH <- function(data, H01, H02, s, u, bl) {
  a <- nrow(H01)
  b <- nrow(H02)
  CH01 <- cumsum(H01[, 3])
  CH02 <- rep(0, a)
  count <- 1
  i = 1
  while (count <= b & i <= a) {
    if (H02[count, 1] < H01[i, 1]) {
      CH02[i] <- CH02[i] + H02[count, 3]
      count <- count + 1
    } else {
      i <- i + 1
    }
  }
  CH02 <- cumsum(CH02)
  
  CIF1 <- 0
  p1a <- nrow(data$Sig) - 1
  b <- as.vector(bl[1:p1a])
  w <- bl[p1a+1]
  for (i in 1:a) {
    if (s < H01[i, 1] & u >= H01[i, 1]) {
      if (i >= 2) {
        CIF1 <- CIF1 + H01[i, 3]*exp(data$X2%*%data$gamma1 + data$alpha1%*%b + data$nu1*w)*
          exp(-CH01[i-1]*exp(data$X2%*%data$gamma1 + data$alpha1%*%b + data$nu1*w)-
                CH02[i-1]*exp(data$X2%*%data$gamma2 + data$alpha2%*%b + data$nu2*w))
      } else {
        CIF1 <- CIF1 + H01[i, 3]*exp(data$X2%*%data$gamma1 + data$alpha1%*%b + data$nu1*w)
      }
      
      
    } else next
  }
  CIF1
}