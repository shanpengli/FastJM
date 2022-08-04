GetVarSurvfit <- function(cdata, ydata, random) {
  
  random.form <- all.vars(random)
  ID <- random.form[length(random.form)]
  n <- nrow(cdata)
  if (length(random.form) == 1) {
    RE <- NULL
    model <- "intercept"
  } else {
    RE <- random.form[-length(random.form)]
    model <- "interslope"
  }
  
  ydim = dim(ydata)
  cdim = dim(cdata)
  
  mdata <- as.data.frame(table(ydata[, ID]))
  colnames(mdata)[1] <- ID
  mdata[, ID] <- as.character(mdata[, ID])
  cdata[, ID] <- as.character(cdata[, ID])
  ydata[, ID] <- as.character(ydata[, ID])
  cmdata <- dplyr::left_join(cdata, mdata, by = ID)
  mdata <- cmdata[, c(1, ncol(cmdata))]
  colnames(mdata) <- c(ID, "ni")
  mdata <- as.data.frame(mdata)
  mdata <- as.vector(mdata$ni)
  mdataS <- rep(0, n) 
  mdataS[1] <- 1
  mdataCum <- cumsum(mdata)
  mdata2 <- mdata - 1
  mdataS[2:n] <- mdataCum[2:n] - mdata2[2:n]
  
  ##random effect covariates
  if (model == "interslope") {
    
    Z <- ydata[, RE]
    Z <- cbind(1, Z)
    Z <- as.matrix(Z)
  } else if (model == "intercept") {
    
    Z <- rep(1, ydim[1])
    Z <- as.data.frame(Z)
    Z <- as.matrix(Z)
    
  } else {
    stop("model should be one of the following options: interslope or intercept.")
  }
  
  ## extract covariates
  X <- ydata[, -c(1:2)]
  X <- as.matrix(cbind(1, X))

  Y <- as.vector(ydata[, 2])
  X2 <- as.matrix(cdata[, -c(1:3)])
  
  survtime <- as.vector(cdata[, 2])
  cmprsk <- as.vector(cdata[, 3])
  
  a <- list(Z = Z, X = X, Y = Y, X2 = X2, survtime = survtime, cmprsk = cmprsk, 
            mdata = mdata, mdataS = mdataS)
  
  return(a)
  
}