logLik <- function(data, b) {
  sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*data$sigma)) + data$CH0*exp(data$X2%*%data$gamma + data$nu%*%b) +
    t(b)%*%solve(data$Sig)%*%b/2
}

logLik.learn <- function(data, b) {
  total <- sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*data$sigma) + 0.5*log(data$sigma)) + 
    data$CH01*exp(data$X2%*%data$gamma1 + data$nu1%*%b) +
    t(b)%*%solve(data$Sig)%*%b/2 + 0.5*log(det(data$Sig))
  if (data$D == 1) {
    total <- total - log(data$HAZ01) - (data$X2%*%data$gamma1 + data$nu1%*%b)
    total
  } else {
    total
  }
}

logLik.JMH <- function(data, bw) {
  p1a <- nrow(data$Sig) - 1
  b <- as.vector(bw[1:p1a])
  w <- bw[p1a+1]
  sigma <- exp(data$W%*%data$tau + w)
  sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*sigma) + 0.5*log(sigma)) + data$CH0*exp(data$X2%*%data$gamma + data$alpha%*%b + data$nu*w) +
    t(bw)%*%solve(data$Sig)%*%bw/2 + 0.5*log(det(data$Sig))
}

logLik.learn.JMH <- function(data, bw) {
  p1a <- nrow(data$Sig) - 1
  b <- as.vector(bw[1:p1a])
  w <- bw[p1a+1]
  sigma <- exp(data$W%*%data$tau + w)
  total <- sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*sigma) + 0.5*log(sigma)) + 
    data$CH01*exp(data$X2%*%data$gamma1 + data$alpha1%*%b + data$nu1*w) +
    t(bw)%*%solve(data$Sig)%*%bw/2 + 0.5*log(det(data$Sig))
  if (data$D == 1) {
    total <- total - log(data$HAZ01) - (data$X2%*%data$gamma1 + data$alpha1%*%b + data$nu1*w)
    total
  } else {
    total
  }
}
