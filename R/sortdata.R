
sortdata <- function(cdata, ydata, ID, surv.formula, long.formula) {
  
  colnames(cdata)[-which(colnames(cdata) == ID)] <- paste0(colnames(cdata)[-which(colnames(cdata) == ID)], ".sur")
  colnames(ydata)[-which(colnames(ydata) == ID)] <- paste0(colnames(ydata)[-which(colnames(ydata) == ID)], ".long")
  
  Tdata <- dplyr::left_join(ydata, cdata, by = ID)
  Truelong <- all.vars(long.formula)
  Truesurvival <- all.vars(surv.formula)
  long <- all.vars(long.formula)
  long <- paste0(long, ".long")
  survival <- all.vars(surv.formula)
  survival <- paste0(survival, ".sur")
  surv <- survival[1]
  
  Tdata <- Tdata[order(-Tdata[, surv], Tdata[, ID]), ]
  
  cdata <- unique(Tdata[, c(ID, survival)])
  ydata <- Tdata[, c(ID, long)]
  colnames(cdata)[-1] <- Truesurvival
  colnames(ydata)[-1] <- Truelong
  
  mdata <- as.data.frame(table(ydata[, ID]))
  colnames(mdata)[1] <- ID
  mdata[, ID] <- as.character(mdata[, ID])
  cdata[, ID] <- as.character(cdata[, ID])
  ydata[, ID] <- as.character(ydata[, ID])
  cmdata <- dplyr::left_join(cdata, mdata, by = ID)
  mdata <- cmdata[, c(1, ncol(cmdata))]
  colnames(mdata) <- c(ID, "ni")
  
  a <- list(ydata, cdata, mdata)
  names(a) <- c("ydata", "cdata", "mdata")
  
  return(a)
  
  
}