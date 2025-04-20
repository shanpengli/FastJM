getmvdummy <- function(long.formula, surv.formula, random, ydata, cdata) {
  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
  cnames <- colnames(cdata)
  ynames <- colnames(ydata)
  
  ID <- all.vars(random)[length(all.vars(random))]
  
  random.var <- all.vars(random)
  
  m <- model.frame(long.formula, ydata)
  ydata2 <- model.matrix(long.formula, m)
  ynames <- colnames(ydata2)
  ydata2 <- data.frame((na.omit(ydata[,c(ID,long)]))[,random.var[length(random.var)]], m[[1]], ydata2[, -1])
  colnames(ydata2) <- c(random.var[length(random.var)], names(m)[1], ynames[-1])
  surv.formula <- surv.formula[-2]
  m <- model.frame(surv.formula, cdata)
  cdata2 <- model.matrix(surv.formula, m)
  cnames <- colnames(cdata2)
  surv.var <- survival[c(1:2)]
  cdata <- cdata[, c(ID, surv.var)]
  cdata2 <- data.frame(cdata, cdata2[, -1])
  colnames(cdata2) <- c(ID, surv.var, cnames[-1])
  
  result <- list(ydata2, cdata2)
  names(result) <- c("ydata", "cdata")
  return(result)
  
}