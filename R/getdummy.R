getdummy <- function(long.formula, surv.formula, random, ydata, cdata) {
  
  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
  cnames <- colnames(cdata)
  ynames <- colnames(ydata)
  ID <- all.vars(random)[length(all.vars(random))]
  
  ##variable check
  if (prod(long %in% ynames) == 0) {
    Fakename <- which(long %in% ynames == FALSE)
    stop(paste0("The variable ", long[Fakename], " not found"))
  }
  if (prod(survival %in% cnames) == 0) {
    Fakename <- which(survival %in% cnames == FALSE)
    stop(paste0("The variable ", survival[Fakename], " not found"))
  }
  if (!(ID %in% ynames)) {
    stop(paste0("ID column ", ID, " not found in the longitudinal dataset!"))
  }
  if (!(ID %in% cnames)) {
    stop(paste0("ID column ", ID, " not found in the survival dataset!"))
  }
  
  
  random.var <- all.vars(random)
  long.var <- long[-1]
  ydata2 <- ydata[,c(random.var[length(random.var)], long)]
  
  long_xnam <- paste(long.var, sep = "")
  long.fixed <- paste(long_xnam, collapse= "+")
  longfmla <- as.formula(paste(" ", long.fixed, sep = "~"))
  test <- as.data.frame(model.matrix(longfmla, data = ydata2[, -1])[, -1])
  ydata2 <- data.frame(ydata2[, 1:2], test)
  long.var <- colnames(ydata2)[-(1:2)]
  long_xnam <- paste(long.var, sep = "")
  long.fixed <- paste(long_xnam, collapse= "+")
  long.formula <- as.formula(paste(colnames(ydata2)[2], long.fixed, sep = "~"))
  
  surv.var <- survival[-c(1:2)]
  cdata2 <- cdata[,c(random.var[length(random.var)], survival)]
  
  surv_xnam <- paste(surv.var, sep = "")
  surv.fixed <- paste(surv_xnam, collapse= "+")
  survfmla <- as.formula(paste(" ", surv.fixed, sep = "~"))
  test <- as.data.frame(model.matrix(survfmla, data = cdata2[, -c(1:3)])[, -1])
  cdata2 <- data.frame(cdata2[, 1:3], test)
  surv.var <- colnames(cdata2)[-(1:3)]
  surv_xnam <- paste(surv.var, sep = "")
  surv.fixed <- paste(surv_xnam, collapse= "+")
  survfmla.out <- paste0("Surv(", survival[1], ", ", survival[2], ")")
  surv.formula <- as.formula(paste(survfmla.out, surv.fixed, sep = "~"))
  
  result <- list(long.formula, ydata2, cdata2, surv.formula)
  names(result) <- c("long.formula", "ydata", "cdata", "surv.formula")
  return(result)
  
}