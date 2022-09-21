GetEmpiricalCIF <- function(data, time, status) {
  library(survival)
  if (!is.data.frame(data))
    stop("data must be a dataframe format.")
  
  data$status1 <- ifelse(data[, status] == 1, 1, 0)
  data$status2 <- ifelse(data[, status] == 2, 1, 0)
  data$time <- data[, time]
  ## Kaplan-Meier estimate S(t) of all-cause failures
  data$status3 <- ifelse(data$status1+data$status2 >= 1, 1, 0)
  fit <- survfit(Surv(time, status3) ~ 1, data = data)
  table <- summary(fit)
  table1 <- data.frame(table$time, table$surv)
  table1 <- rbind(c(0, 1), table1)
  table1$table.surv <- c(1, table1$table.surv[1:(nrow(table1)-1)])
  colnames(table1) <- c("time", "survProb")
  
  ## NA estimate of H1(t)
  fit <- survfit(Surv(time, status1) ~ 1, data = data)
  tableH1 <- summary(fit)
  h1 <- tableH1$n.event/tableH1$n.risk
  time1 <- tableH1$time
  H1 <- data.frame(time1, h1)
  ## CIF estimate of risk 1
  SurvProb <- vector()
  CIF1 <- vector()
  CIF <- 0
  count <- 1
  for (i in 1:nrow(table1)) {
    if (count > length(time1)) {
      break
    } else if (table1[i, 1] != time1[count]) {
      next
    } else {
      CIF <- CIF + table1[i, 2]*h1[count]
      CIF1[count] <- CIF
      SurvProb[count] <- table1[i, 2]
      count <- count + 1
    }
  }
  H1$SurvProb <- SurvProb
  H1$CIF1 <- CIF1
  
  ## NA estimate of H2(t)
  fit <- survfit(Surv(time, status2) ~ 1, data = data)
  tableH2 <- summary(fit)
  h2 <- tableH2$n.event/tableH2$n.risk
  time2 <- tableH2$time
  H2 <- data.frame(time2, h2)
  ## CIF estimate of risk 2
  SurvProb <- vector()
  CIF2 <- vector()
  CIF <- 0
  count <- 1
  for (i in 1:nrow(table1)) {
    if (count > length(time2)) {
      break
    } else if (table1[i, 1] != time2[count]) {
      next
    } else {
      CIF <- CIF + table1[i, 2]*h2[count]
      CIF2[count] <- CIF
      SurvProb[count] <- table1[i, 2]
      count <- count + 1
    }
  }
  H2$SurvProb <- SurvProb
  H2$CIF2 <- CIF2
  
  return(list(H1 = H1, H2 = H2))
}