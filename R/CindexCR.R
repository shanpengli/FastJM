CindexCR <- function (time, status, predicted, Cause_int = 1, Time = NULL) {
  if (any(is.na(time))) {
    stop("The input vector cannot have NA")
  }
  if (any(is.na(status))) {
    stop("The input vector cannot have NA")
  }
  if (any(!(status %in% c(0, 1, 2)))) {
    stop("The status must be 0 or 1 or 2")
  }
  if (any(is.na(predicted))) {
    stop("The input vector cannot have NA")
  }
  if (!(Cause_int %in% status)) {
    stop("Invalid input of Cause_int")
  }
  if (min(time) <= 0) {
    stop("Survival time must be positive")
  }
  
  KM <- summary(survfit(Surv(time, status==0)~1))
  df_KM <- data.frame(time = c(0, KM$time), surv = c(1, KM$surv))
  
  Time_survival <- time
  Censoring <- ifelse(status == 0, 0, 1)
  Cause <- ifelse(status == 2, 2, 1)
  Prediction <- log(predicted)
  if(is.null(Time)) Time <- max(Time_survival) + 1
  n <- length(Prediction)
  A <- matrix(0, nrow = n, ncol = n)
  B <- matrix(0, nrow = n, ncol = n)
  Q <- matrix(0, nrow = n, ncol = n)
  N_t <- matrix(0, nrow = n, ncol = n)
  Num_mat <- matrix(0, nrow = n, ncol = n)
  Den_mat <- matrix(0, nrow = n, ncol = n)
  Num <- 0
  Den <- 0
  
  w1 <- rep(0, n)
  w2 <- rep(0, n)
  for(i in 1:n){
    w1[i] <- min(df_KM$surv[df_KM$time < Time_survival[i]]) 
    w2[i] <- min(df_KM$surv[df_KM$time <= Time_survival[i]])
  }
  w3 <- w1 %*% t(w1)
  w1 <- w1 * w2
  
  for (i in 1:n) {
    A[i, which(Time_survival[i] < Time_survival)] <- 1
    A[i, ] <- A[i, ] * w1[i]
    B[i, intersect(intersect(which((Time_survival[i] >= 
                                      Time_survival)), which(Cause != Cause_int)), 
                   which(Censoring == 1))] <- 1
    B[i, ] <- B[i, ] * w3[i, ]
    Q[i, which(Prediction[i] > Prediction)] <- 1
  }
  for (i in 1:n) {
    if (Time_survival[i] <= Time && Cause[i] == Cause_int && 
        Censoring[i] == 1) {
      N_t[i, ] <- 1
    }
  }
  Num_mat <- (A + B) * Q * N_t
  Den_mat <- (A + B) * N_t
  Num <- sum(Num_mat)
  Den <- sum(Den_mat)
  return(Num/Den)
}