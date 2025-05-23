sortmvdata <- function(cdata, ydata, ID, surv.formula, long.formula) {
  
  colnames(cdata)[-which(colnames(cdata) == ID)] <- paste0(colnames(cdata)[-which(colnames(cdata) == ID)], ".sur") # add ".sur" to column names in cdata except for ID
  colnames(ydata)[-which(colnames(ydata) == ID)] <- paste0(colnames(ydata)[-which(colnames(ydata) == ID)], ".long") # add ".long" to column names in ydata except for ID
  
  Tdata <- dplyr::left_join(ydata, cdata, by = ID) # merge ydata and cdata by ID
  
  Truelong <- long <- c()
  
  numBio <- length(long.formula)
  
  for(g in 1:numBio){
    Truelong <- c(Truelong,all.vars(long.formula[[g]])) # get variable names from the long formula
    long <- c(long, all.vars(long.formula[[g]])) # prepare column names for long data
  }
  
  Truelong <- unique(Truelong)
  long <- unique(long)
  long <- paste0(long, ".long")
  
  Truesurvival <- all.vars(surv.formula) # get variable names from the survival formula
  
  survival <- all.vars(surv.formula) # prepare column names for survival data
  survival <- paste0(survival, ".sur")
  surv <- survival[1] # choose the first survival variable
  
  Tdata <- Tdata[order(-Tdata[, surv], Tdata[, ID]), ] # sort merged data by survival value (descending) and then by ID
  
  cdata <- unique(Tdata[, c(ID, survival)]) # create a unique dataset for cdata with only ID and survival columns
  
  ## this line problematic
  ydata <- Tdata[, c(ID, long)] # create a dataset for ydata with ID and long columns
  colnames(cdata)[-1] <- Truesurvival # rename cdata columns to original survival variable names
  colnames(ydata)[-1] <- Truelong # rename ydata columns to original long variable names
  
  
  cdata[, ID] <- as.character(cdata[, ID])
  
  ydataList <- mdataList <- vector("list",numBio)
  for(g in 1:numBio){
    ydatatemp <- na.omit(ydata[, c(ID, all.vars(long.formula[[g]]))])
    mdatatemp <- as.data.frame(table(ydatatemp[, ID])) # create a count of occurrences of each ID in ydata
    colnames(mdatatemp)[1] <- ID
    mdatatemp[,ID] <- as.character(mdatatemp[, ID])
    ydatatemp[, ID] <- as.character(ydatatemp[, ID])
    
    cmdata <- dplyr::left_join(cdata, mdatatemp, by = ID)
    mdatatemp <- cmdata[, c(1, ncol(cmdata))]
    colnames(mdatatemp) <- c(ID, "ni") # rename count column to "ni"
    
    
    ydataList[[g]] <- ydatatemp
    mdataList[[g]] <- mdatatemp
    
    
  }
  
  
  # mdata <- as.data.frame(table(ydata[, ID])) # create a count of occurrences of each ID in ydata
  # colnames(mdata)[1] <- ID
  # mdata[, ID] <- as.character(mdata[, ID])
  # cdata[, ID] <- as.character(cdata[, ID])
  # ydata[, ID] <- as.character(ydata[, ID])
  # cmdata <- dplyr::left_join(cdata, mdata, by = ID) # merge cdata and count data by ID
  # mdata <- cmdata[, c(1, ncol(cmdata))]
  # colnames(mdata) <- c(ID, "ni") # rename count column to "ni"
  
  a <- list(ydataList, cdata, mdataList) # return a list with ydata, cdata, and mdata
  names(a) <- c("ydata", "cdata", "mdata")
  
  return(a)
  
}
