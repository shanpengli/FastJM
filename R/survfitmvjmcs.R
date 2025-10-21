##' @title Prediction in Joint Models
##' @name survfitmvjmcs
##' @aliases survfitmvjmcs
##' @description This function computes the conditional probability of 
##' surviving later times than the last observed time for which a longitudinal 
##' measurement was available.
##' @param object an object inheriting from class \code{mvjmcs}.
##' @param seed a random seed number to proceed Monte Carlo simulation. Default is 100.
##' @param ynewdata a data frame that contains the longitudinal and covariate information for the subjects 
##' for which prediction of survival probabilities is required.
##' @param cnewdata a data frame that contains the survival and covariate information for the subjects 
##' for which prediction of survival probabilities is required.
##' @param u a numeric vector of times for which prediction survival probabilities are to be computed.
##' @param Last.time a numeric vector or character string. This specifies the known time at which each of 
##' the subjects in cnewdata was known to be alive. If NULL, then this is automatically taken as the 
##' survival time of each subject. If a numeric vector, then it is assumed to be greater than or equals to the 
##' last available longitudinal time point for each subject. If a character string, then it should be 
##' a variable in cnewdata.
##' @param obs.time a character string of specifying a longitudinal time variable in ynewdata.
##' @param LOCF a logical value to indicate whether the last-observation-carried-forward approach applies to prediction. 
##' If \code{TRUE}, then \code{LOCFcovariate} and \code{clongdata} must be specified to indicate 
##' which time-dependent survival covariates are included for dynamic prediction. Default is FALSE.
##' @param LOCFcovariate a vector of string with time-dependent survival covariates if \code{LOCF = TRUE}. Default is NULL.
##' @param clongdata a long format data frame where time-dependent survival covariates are incorporated. Default is NULL.
##' @param ... further arguments passed to or from other methods. 
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{mvjmcs}}
##' @export
survfitmvjmcs <- function(object, seed = 100, ynewdata = NULL, cnewdata = NULL, 
                        u = NULL, Last.time = NULL, obs.time = NULL, 
                        LOCF = FALSE, LOCFcovariate = NULL, clongdata = NULL, ...) {
  
  
  # x2 --> 2
  # nu --> alpha
  if (!inherits(object, "mvjmcs"))
    stop("Use only with 'mvjmcs' objects.\n")
  if (is.null(ynewdata))
    stop("New longitudinal data for dynamic prediction is needed.")
  if (is.null(cnewdata))
    stop("New longitudinal data for dynamic prediction is needed.")
  if (is.null(u)) 
    stop("Please specify the future time for dynamic prediction.")   
  # if (!method %in% c("Laplace", "GH"))
  #   stop("Please specify a method for probability approximation: Laplace or GH.")
  if (!is.vector(u)) 
    stop("u must be vector typed.")
  # if (is.null(quadpoint)) {
  #   quadpoint <- object$quadpoint
  # }
  if (is.null(obs.time)) {
    stop("Please specify a vector that represents the time variable from ydatanew.")
  } else {
    if (!obs.time %in% colnames(ynewdata)) {
      stop(paste0(obs.time, " is not found in ynewdata."))
    }
  }
  
  
  # check number of Bio
  pRE <- c()
  if(!is.list(object$random)){
    bvar <- all.vars(object$random)
    pRE[1] <- length(all.vars(object$random))
    ID <- bvar[length(bvar)]
    numBio = 1
    
  } else{
    numBio = length(object$random)
    bvar <- vector("list", numBio)
    
    for(g in 1:numBio){
      bvar[[g]] <- all.vars(object$random[[g]])
      pRE[g] <- length(all.vars(object$random[[g]]))
    }
      ID <- bvar[[g]][length(bvar[[g]])]
  }
  
  if (!(ID %in% colnames(ynewdata)))
    stop(paste("The ID variable", ID, "is not found in ynewdata."))
  if (!(ID %in% colnames(cnewdata)))
    stop(paste("The ID variable", ID, "is not found in cnewdata."))
  
  if (LOCF) {
    if (is.null(clongdata)) {
      stop("Please provide a long format data frame that includes all longitudinal covariates for dynamic prediction.\n")
    } else {
      survival <- all.vars(object$SurvivalSubmodel) # no Survival submodel
      if (prod(LOCFcovariate %in% survival) == 0)
        stop("Some covariates are not trained in the joint model. Please reconsider the covariates for prediction.\n")
      if (!obs.time %in% colnames(clongdata))
        stop(paste0(obs.time, " is not found in clongdata."))
      if (!(ID %in% colnames(clongdata)))
        stop(paste("The ID variable", ID, "is not found in clongdata."))
      clongID <- unique(clongdata[, ID])
      cID <- cnewdata[, ID]
      if (prod(clongID == cID) == 0) {
        stop("The order of subjects in clongdata doesn't match with cdata.")
      }
    }
  }
  
  # if(is.null(stime)){
  #   stime = 0
  # }

  if (LOCF) {
    clongdata <- clongdata[clongdata[, obs.time] <= Last.time, ]
    clongdataLOCF <- clongdata %>%
      dplyr::group_by(dplyr::across(ID)) %>%
      dplyr::filter(row_number() == n())
    cnewdata[, LOCFcovariate] <- clongdataLOCF[, LOCFcovariate]
  }
  
  # get only relevant variables from new dataset
  ynewdata <- ynewdata[, colnames(object$ydata)]
  cnewdata <- cnewdata[, colnames(object$cdata)]
  
  # merge datasets together
  ydata2 <- rbind(object$ydata, ynewdata)
  cdata2 <- rbind(object$cdata, cnewdata)
  
  # if statement for checking numbio matchess
  if(is.list(object$LongitudinalSubmodel)){
    numBio <- length(object$LongitudinalSubmodel)
  }else{
    numBio <- 1
  }
  
  
  ### CHECK IF NOT LIST/IF BIO = 1
  
  ydata2List <- ynewdata2 <- ynewdatasplit <- vector("list", numBio)
  ny <- Ny <- c()
  
  Yvar <-  vector("list", numBio)
  
  # Yvar <- colnames(ydata2)[-1]
  # Cvar <- colnames(cdata2)[-1]
  
  cnewdata <- cnewdata[, colnames(object$cdata)]
  cdata2 <- rbind(object$cdata, cnewdata)
  
  nc <- nrow(cnewdata)
  
  if(numBio == 1){
    
    bvar1 <- vector("list", 1)
    if (length(bvar) > 1) {
      bvar1[[1]] <- bvar[-length(bvar)]
    } else {
      bvar1[[1]] <- character(0)
    }
    
    getdum <- getdummy(long.formula = object$LongitudinalSubmodel,
                       surv.formula = object$SurvivalSubmodel,
                       random = object$random, ydata = ydata2, cdata = cdata2)
    
    g = 1
  
    ydata2List[[g]] <- getdum$ydata
    ydata2temp <- getdum$ydata
    cdata2temp <- getdum$cdata
    ynewdatasplit[[g]] <- na.omit(ynewdata[,all.vars(object$LongitudinalSubmodel)])
    
    
    Yvar[[g]] <- colnames(ydata2temp)[-1]
    Cvar <- colnames(cdata2temp)[-1]
    # bvar <- all.vars(object$random) # CHANGE THIS PART TO ACCOUNT FOR MULTIPLE BIOMARKERS
    
    ny[g] <- nrow(ynewdatasplit[[g]])
    Ny[g] <- nrow(ydata2temp)
    
    ynewdata2[[g]] <- ydata2temp[c((Ny[g]-ny[g]+1):Ny[g]), ]
    
    }else{
  
      bvar1 <- vector("list", numBio)
      for(g in 1:numBio){

    # object$random needs to be list too
    

      getdum <- getmvdummy(long.formula = object$LongitudinalSubmodel[[g]],
                           surv.formula = object$SurvivalSubmodel,
                           random = object$random[[g]], ydata = ydata2, cdata = cdata2)
    
      ydata2List[[g]] <- getdum$ydata
      ydata2temp <- getdum$ydata
      cdata2temp <- getdum$cdata
      ynewdatasplit[[g]] <- na.omit(ynewdata[,all.vars(object$LongitudinalSubmodel[[g]])])
      
      Yvar[[g]] <- colnames(ydata2temp)[-1]
      Cvar <- colnames(cdata2temp)[-1]

      if (length(bvar[[g]]) > 1) {
        # RE covariates (drop id)
        bvar1[[g]] <- bvar[[g]][-(length(bvar[[g]]))]
      } else {
        # no RE covariates
        bvar1[[g]] <- character(0)
      }
        
      ny[g] <- nrow(ynewdatasplit[[g]])
      Ny[g] <- nrow(ydata2temp)
        
    
      ynewdata2[[g]] <- ydata2temp[c((Ny[g]-ny[g]+1):Ny[g]), ]
    
    }
  }
  
  Nc <- nrow(cdata2temp)
  cnewdata2 <- cdata2temp[c((Nc-nc+1):Nc), ]
  
  survival <- all.vars(object$SurvivalSubmodel)
  cmprsk <- as.vector(cnewdata2[, survival[2]])
  
  nsig <- nrow(object$Sig)
  
  cID <- cnewdata2[, ID]
  yID <- vector("list", numBio)
  for(g in 1:numBio){
    yID[[g]] <- unique((ynewdata2[[g]])[, ID])
  }
  commonyID <- Reduce(intersect, yID)
  
  cnewdata2 <- cnewdata2 %>%
    dplyr::filter(.data[[ID]] %in% commonyID)
  
  for(g in 1:numBio){
    ynewdata2[[g]] <- ynewdata2[[g]] %>%
      dplyr::filter(.data[[ID]] %in% commonyID)
  }
  N.ID <- length(commonyID)
  
  if (!is.null(Last.time)) {
    if (is.character(Last.time)) {
      if (Last.time %in% colnames(cnewdata2)) {
        Last.time <- cnewdata2[, Last.time]
      } else {
        stop(paste(Last.time, "is not found in cnewdata."))
      }
    }
    if (is.numeric(Last.time)) {
      if (length(Last.time) == 1) {
        Last.time <- rep(Last.time, nrow(cnewdata2))
      } else if (length(Last.time) < nrow(cnewdata2)) {
        stop("The last.time vector does not match cnewdata.")
      } else {
      }
    }
  } else {
    Last.time <- cnewdata2[, Cvar[1]]
  }
 
  Pred <- list()
  y.obs <- list()
  CompetingRisk <- object$CompetingRisk
  
  if (CompetingRisk) {
    # betaFlat <- object$beta
    # betaList <- vector("list",numBio)
    # index <- 1
    # for(g in 1:numBio){
    #   p <- length(all.vars(object$LongitudinalSubmodel[[g]]))
    #   betaList[[g]] <- betaFlat[index:(index + p-1)]
    #   index <- index + g
    # }
    betaList <- object$betaList
    sigma <- object$sigma
    gamma1 <- object$gamma1
    gamma2 <- object$gamma2
    
    alpha1flat <- object$alpha1
    alpha2flat <- object$alpha2
    alpha1 <- alpha2 <- vector("list", numBio)
    index = 1
    for(g in 1:numBio){
      alpha1[[g]] <- alpha1flat[index:(index+pRE[g]-1)] 
      alpha2[[g]] <- alpha2flat[index:(index+pRE[g]-1)] 
      index = index + pRE[g]
    }
    
    H01 <- object$H01
    H02 <- object$H02
    Sig <- object$Sig
    
    Predraw1 <- matrix(0, nrow = nrow(cnewdata2), ncol = length(u))
    Predraw2 <- matrix(0, nrow = nrow(cnewdata2), ncol = length(u))
    lengthu <- length(u)
   
    subNDy <- YList <- XList <- y.obs <-  ZList<- vector("list", numBio)
    
    # n <- nrow(cnewdata2)
    

    for (j in 1:N.ID) {
      for(g in 1:numBio){
        subNDy[[g]] <- ynewdata2[[g]][ynewdata2[[g]][, ID] == commonyID[j], ]
        y.obs[[g]][[j]] <- data.frame(ynewdata[ynewdata[, ID] == commonyID[j], c(obs.time, Yvar[[g]][1])])
        YList[[g]] <- subNDy[[g]][, Yvar[[g]][1]]
        XList[[g]] <- as.matrix(data.frame(1, subNDy[[g]][, Yvar[[g]][-1]]))
        if (pRE[g] == 1) {
          ZList[[g]] <- matrix(1, ncol = 1, nrow = length(YList[[g]]))
        } else {
          ZList[[g]] <- as.matrix(data.frame(1, subNDy[[g]][, bvar1[[g]], drop = FALSE]))
        }
        
      }

      subNDc <- cnewdata2[cnewdata2[, ID] == commonyID[j], ]
    
      stime <-  as.numeric(Last.time[j])
      CH01 <- CH(H01, stime)
      CH02 <- CH(H02, stime)

      W <- as.matrix(subNDc[1, Cvar[3:length(Cvar)]])

  
    
    data <- list(YList, XList, ZList, W = W, CH01 = CH01, CH02 = CH02, betaList, gamma1, gamma2, alpha1, alpha2, sigma, Sig)
    names(data) <- c("YList", "XList", "ZList", "W", "CH01", "CH02", "beta",
                     "gamma1", "gamma2", "alpha1", "alpha2",  "sigma", "Sig")
    
    opt <- optim(rep(0, nsig), logLikCRmv, data = data, method = "BFGS", hessian = TRUE)
    meanb <- opt$par
    Poscov <- solve(opt$hessian)
      
      
        pREvec <- c()
        
        for(g in 1:numBio){
          pREvec[g] <- ncol(ZList[[g]])
        }
        
        for (jj in 1:lengthu) {
          ## calculate the CIF
          CIF1 <- CIF1mv.CR(data, H01, H02, stime, u[jj], meanb, numBio, pREvec)
          P1us <- Pkmv.us(CIF1, data, meanb, numBio = numBio, pREvec)
          Predraw1[j, jj] <- P1us
          
          CIF2 <- CIF2mv.CR(data, H01, H02, stime, u[jj], meanb, numBio, pREvec)
          P2us <- Pkmv.us(CIF2, data, meanb, numBio = numBio, pREvec)
          Predraw2[j, jj] <- P2us
        }
        
    } 
    for (jj in 1:N.ID) {
      Pred[[jj]] <- data.frame(u, Predraw1[jj, ], Predraw2[jj, ])
      colnames(Pred[[jj]]) <- c("times", "CIF1", "CIF2")
    }
    } else{
    
      betaList <- object$betaList
      sigma <- object$sigma
      gamma1 <- object$gamma1
      
      alpha1flat <- object$alpha1
      alpha1 <- vector("list", numBio)
      index = 1
      for(g in 1:numBio){
        alpha1[[g]] <- alpha1flat[index:(index+pRE[g]-1)] 
        index = index + pRE[g]
      }
      
      H01 <- object$H01
      Sig <- object$Sig
      
      Predraw1 <- matrix(0, nrow = nrow(cnewdata2), ncol = length(u))
      lengthu <- length(u)
      
      subNDy <- YList <- XList <- y.obs <-  ZList <- vector("list", numBio)
      
      # n <- nrow(cnewdata2)
      
      
      for (j in 1:N.ID) {
        for(g in 1:numBio){
          subNDy[[g]] <- ynewdata2[[g]][ynewdata2[[g]][, ID] == commonyID[j], ]
          y.obs[[g]][[j]] <- data.frame(ynewdata[ynewdata[, ID] == commonyID[j], c(obs.time, Yvar[[g]][1])])
          YList[[g]] <- subNDy[[g]][, Yvar[[g]][1]]
          XList[[g]] <- as.matrix(data.frame(1, subNDy[[g]][, Yvar[[g]][-1]]))
          if (pRE[g] == 1) {
            ZList[[g]] <- matrix(1, ncol = 1, nrow = length(YList[[g]]))
          } else {
            ZList[[g]] <- as.matrix(data.frame(1, subNDy[[g]][, bvar1[[g]], drop = FALSE]))
          }
          
        }
        
        subNDc <- cnewdata2[cnewdata2[, ID] == commonyID[j], ]
        
        stime <-  as.numeric(Last.time[j])
        CH01 <- CH(H01, stime)
        
        W <- as.matrix(subNDc[1, Cvar[3:length(Cvar)]])
        
        
        data <- list(YList, XList, ZList, W = W, CH01 = CH01,  betaList, gamma1,  alpha1,  sigma, Sig)
        names(data) <- c("YList", "XList", "ZList", "W", "CH01",  "beta",
                         "gamma1","alpha1",   "sigma", "Sig")
        
        opt <- optim(rep(0, nsig), logLikSFmv, data = data, method = "BFGS", hessian = TRUE)
        meanb <- opt$par
        Poscov <- solve(opt$hessian)
        
        pREvec <- c()
        
        for(g in 1:numBio){
          pREvec[g] <- ncol(ZList[[g]])
        }
        
        for (jj in 1:lengthu) {
          ## calculate the CIF
          CIF1 <- CIF1mv.SF(data, H01, stime, u[jj], meanb, numBio, pREvec)
          P1us <- Pkmv.us_SF(CIF1, data, meanb, numBio = numBio, pREvec)
          Predraw1[j, jj] <- 1 - P1us
        }
        
      } 
      
      for (jj in 1:N.ID) {
        Pred[[jj]] <- data.frame(u, Predraw1[jj, ])
        colnames(Pred[[jj]]) <- c("times", "PredSurv")
      }
  }
  
  for(g in 1:numBio){
    names(y.obs[[g]]) <- commonyID
  }
  
  names(Pred) <- commonyID
  Last.time <- data.frame(commonyID, Last.time)
  colnames(Last.time) <- c(ID, "Last.time")
  sum <- list(Pred = Pred, Last.time = Last.time, y.obs = y.obs,
              CompetingRisk = CompetingRisk)
  class(sum) <- "survfitmvjmcs"
  sum
}

