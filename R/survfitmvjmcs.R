##' @title Dynamic predictions from mvjmcs models
##' @name survfitmvjmcs
##' @aliases survfitmvjmcs
##' @description This function computes the conditional probability of 
##' surviving later times than the last observed time for which a longitudinal 
##' measurement was available.
##' 
##' @details Internal implementation used by \code{\link{survfitJM}}.
##' Call \code{survfitJM()} for dynamic predictions.
##'
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
##' @keywords internal
survfitmvjmcs <- function(object, seed = 100, ynewdata = NULL, cnewdata = NULL, 
                          u = NULL, Last.time = NULL, obs.time = NULL, 
                          LOCF = FALSE, LOCFcovariate = NULL, clongdata = NULL, ...) {
  
  # x2 --> 2
  # nu --> alpha
  
  # validate inputs
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
  
  # check longitudinal time variable
  if (is.null(obs.time)) {
    stop("Please specify a vector that represents the time variable from ydatanew.")
  } else {
    if (!obs.time %in% colnames(ynewdata)) {
      stop(paste0(obs.time, " is not found in ynewdata."))
    }
  }
  
  
  # check number of biomarkers and RE per biomarker
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
  
  # check subject ID variable in new datasets
  if (!(ID %in% colnames(ynewdata)))
    stop(paste("The ID variable", ID, "is not found in ynewdata."))
  if (!(ID %in% colnames(cnewdata)))
    stop(paste("The ID variable", ID, "is not found in cnewdata."))
  
  survival <- all.vars(object$SurvivalSubmodel)
  
  
  latAsso <- object$latAsso
  
  if(latAsso == "sre"){
    landmark  <- FALSE
  }else{
    landmark <- TRUE
  }
  
  if (!latAsso %in% c("sre", "present", "presentlp")) {
    stop("latAsso must be one of 'sre', 'present', or 'presentlp'.")
  }
  
  
  
  # determine prediction time for each subject
  if (!is.null(Last.time)) {
    if (is.character(Last.time)) {
      if (Last.time %in% colnames(cnewdata)) {
        Last.time <- cnewdata[, Last.time]
      } else {
        stop(paste(Last.time, "is not found in cnewdata."))
      }
    }
    
    if (is.numeric(Last.time)) {
      if (length(Last.time) == 1) {
        Last.time <- rep(Last.time, nrow(cnewdata))
      } else if (length(Last.time) != nrow(cnewdata)) {
        stop("The Last.time vector does not match cnewdata.")
      }
    }
  } else {
    
    # default to observed survival time
    Last.time <- cnewdata[, survival[1]]
  }
  
  
  
  
  
  names(Last.time) <- as.character(cnewdata[, ID])
  
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
    }
  }
  
  
  
  
  
  if (LOCF) {
    # subset observations up to each subject's prediction time
    clongdata <- clongdata[which(clongdata[, obs.time] <= Last.time[as.character(clongdata[, ID])]), ]
    
    # select the most recent observation for each subject
    clongdataLOCF <- clongdata %>%
      dplyr::group_by(dplyr::across(ID)) %>%
      dplyr::arrange(.data[[obs.time]]) %>%
      dplyr::filter(dplyr::row_number() == dplyr::n()) %>%
      dplyr::ungroup()
    
    # match subjects by ID
    index <- match(cnewdata[, ID], clongdataLOCF[[ID]])
    
    if (anyNA(index)) {
      stop("No LOCF observations available by Last.time.")
    }
    
    cnewdata[, LOCFcovariate] <- clongdataLOCF[index, LOCFcovariate]
  }
  
  # retain relevant variables
  ynewdata <- ynewdata[, colnames(object$ydata)]
  cnewdata <- cnewdata[, colnames(object$cdata)]
  
  # merge training and pred datasets together for dummy
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
  
  cnewdata <- cnewdata[, colnames(object$cdata)]
  cdata2 <- rbind(object$cdata, cnewdata)
  
  nc <- nrow(cnewdata)
  
  
  # ~~~~~~~~~~~~~~~~
  # Single Biomarker
  # ~~~~~~~~~~~~~~~~
  
  if(numBio == 1){
    
    bvar1 <- vector("list", 1)
    if (length(bvar) > 1) {
      bvar1[[1]] <- bvar[-length(bvar)]
    } else {
      bvar1[[1]] <- character(0)
    }
    
    # dummy encoding
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
    
    # count prediction obs and combined obs
    ny[g] <- nrow(ynewdatasplit[[g]])
    Ny[g] <- nrow(ydata2temp)
    
    ynewdata2[[g]] <- ydata2temp[c((Ny[g]-ny[g]+1):Ny[g]), ]
    
  }else{
    
    bvar1 <- vector("list", numBio)
    for(g in 1:numBio){
      
      # object$random needs to be list too
      
      # dummy encoding
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
  
  
  commonyID <- intersect(Reduce(intersect, yID), cID)
  
  # retain eligible subjects
  cnewdata2 <- cnewdata2 %>% dplyr::filter(.data[[ID]] %in% commonyID)
  for(g in 1:numBio){
    ynewdata2[[g]] <- ynewdata2[[g]] %>% dplyr::filter(.data[[ID]] %in% commonyID)
  }
  N.ID <- length(commonyID)
  
  if (N.ID == 0) {
    stop("No common subjects found in the prediction datasets.")
  }
  
  # if (landmark) {
  #   lm.time <- object$s
  # }
  
  Last.time <- unname(Last.time[as.character(commonyID)])
  
  
  Pred <- list()
  y.obs <- list()
  CompetingRisk <- object$CompetingRisk
  
  # ~~~~~~~~~~~~~~~
  # Competing Risks
  # ~~~~~~~~~~~~~~~
  
  if (CompetingRisk) {
    
    # initialization
    betaList <- object$betaList
    sigma <- object$sigma
    gamma1 <- object$gamma1
    gamma2 <- object$gamma2
    
    alpha1flat <- object$alpha1
    alpha2flat <- object$alpha2
    
    if (latAsso == "sre") {
      
      alpha1 <- alpha2 <- vector("list", numBio)
      index <- 1
      
      for (g in 1:numBio) {
        alpha1[[g]] <- alpha1flat[index:(index + pRE[g] - 1)]
        alpha2[[g]] <- alpha2flat[index:(index + pRE[g] - 1)]
        index <- index + pRE[g]
      }
      
    } else {
      
      # present/presentlp: one association parameter per biomarker
      alpha1 <- alpha1flat
      alpha2 <- alpha2flat
      
    }
    
    H01 <- object$H01
    H02 <- object$H02
    Sig <- object$Sig
    
    Predraw1 <- matrix(0, nrow = nrow(cnewdata2), ncol = length(u))
    Predraw2 <- matrix(0, nrow = nrow(cnewdata2), ncol = length(u))
    lengthu <- length(u)
    
    subNDy <- YList <- XList <- y.obs <-  ZList<- vector("list", numBio)
    
    # n <- nrow(cnewdata2)
    
    # perform dynamic prediction per subject
    for (j in 1:N.ID) {
      
      # time s which subject last known to be event-free!!
      stime <-  as.numeric(Last.time[j])
      
      for(g in 1:numBio){
        
        # if(landmark){
        #   # use longitudinal history available by landmark
        #   subNDy[[g]] <- ynewdata2[[g]][ynewdata2[[g]][, ID] == commonyID[j] & ynewdata2[[g]][, obs.time] <= stime,,]
        #   y.obs[[g]][[j]] <- data.frame(ynewdata[ynewdata[, ID] == commonyID[j] & ynewdata[, obs.time] <= stime,
        #                                          c(obs.time, Yvar[[g]][1])])
        #   } else{
        #     subNDy[[g]] <- ynewdata2[[g]][ynewdata2[[g]][, ID] == commonyID[j], ]
        #     y.obs[[g]][[j]] <- data.frame(ynewdata[ynewdata[, ID] == commonyID[j], c(obs.time, Yvar[[g]][1])])
        # }
        
        # retain longitudinal history available at or before time s
        subNDy[[g]] <- ynewdata2[[g]][ynewdata2[[g]][, ID] == commonyID[j] &
                                        ynewdata2[[g]][, obs.time] <= stime, , drop = FALSE]
        
        # store obs biomarker measurements up to time s
        y.obs[[g]][[j]] <- data.frame(ynewdata[ynewdata[, ID] == commonyID[j] & 
                                                 ynewdata[, obs.time] <= stime,
                                               c(obs.time, Yvar[[g]][1])])
        
        # order longitudinal observations by time
        subNDy[[g]] <- subNDy[[g]][order(subNDy[[g]][, obs.time]), ,drop = FALSE]
        
        
        if (nrow(subNDy[[g]]) == 0) {
          stop("No longitudinal observations available by Last.time.")
        }
        
        # matrix construction for posterior est.
        YList[[g]] <- subNDy[[g]][, Yvar[[g]][1]]
        XList[[g]] <- as.matrix(data.frame(1, subNDy[[g]][, Yvar[[g]][-1]]))
        if (pRE[g] == 1) {
          # intercept only
          ZList[[g]] <- matrix(1, ncol = 1, nrow = length(YList[[g]]))
        } else {
          # intercept + RE
          ZList[[g]] <- as.matrix(data.frame(1, subNDy[[g]][, bvar1[[g]], drop = FALSE]))
        }
        
      }
      
      subNDc <- cnewdata2[cnewdata2[, ID] == commonyID[j], ]
      
      # matrix construction for non-SRE; initalized at landmark time
      
      Xs_i <- NULL
      Zs_i <- NULL
      
      if (landmark && latAsso %in% c("present", "presentlp")) {
        
        Xs_i <- vector("list", numBio)
        Zs_i <- vector("list", numBio)
        
        for (g in 1:numBio) {
          
          # fixed-effects design at prediction time
          if (latAsso == "present") {
            
            Xs_i[[g]] <- XList[[g]][nrow(XList[[g]]), , drop = FALSE]
            
            time_fe <- which(Yvar[[g]][-1] == obs.time) + 1
            
            if (length(time_fe) > 0) {
              Xs_i[[g]][1, time_fe] <- stime
            }
          }
          
          # random-effects design at prediction time
          if (pRE[g] == 1) {
            
            Zs_i[[g]] <- matrix(1, nrow = 1, ncol = 1)
            
          } else {
            
            Zs_i[[g]] <- ZList[[g]][nrow(ZList[[g]]), , drop = FALSE]
            
            time_re <- which(bvar1[[g]] == obs.time) + 1
            
            if (length(time_re) > 0) {
              Zs_i[[g]][1, time_re] <- stime
            }
          }
        }
      }
      
      CH01 <- CH(H01, stime)
      CH02 <- CH(H02, stime)
      
      W <- as.matrix(subNDc[1, Cvar[3:length(Cvar)]])
      
      data <- list(
        YList = YList,
        XList = XList,
        ZList = ZList,
        W = W,
        CH01 = CH01,
        CH02 = CH02,
        beta = betaList,
        gamma1 = gamma1,
        gamma2 = gamma2,
        alpha1 = alpha1,
        alpha2 = alpha2,
        sigma = sigma,
        Sig = Sig,
        latAsso = latAsso,
        landmark = landmark,
        H01 = H01,
        H02 = H02,
        stime = stime,
        Xs_i = Xs_i,
        Zs_i = Zs_i
      )
      
      
      # use the last guess from fitting as initial
      opt <- optim(rep(0, nsig), logLikCRmv, data = data, method = "BFGS")
      meanb <- opt$par
      
      pREvec <- c()
      for(g in 1:numBio){
        pREvec[g] <- ncol(ZList[[g]])
      }
      
      # calculate conditional CIFs at time u
      for (jj in 1:lengthu) {
        if (u[jj] > stime) {
          ## calculate the CIF
          CIF1 <- CIF1mv.CR(data, H01, H02, stime, u[jj], meanb, numBio, pREvec)
          P1us <- Pkmv.us(CIF1, data, meanb, numBio = numBio, pREvec)
          Predraw1[j, jj] <- P1us
          
          CIF2 <- CIF2mv.CR(data, H01, H02, stime, u[jj], meanb, numBio, pREvec)
          P2us <- Pkmv.us(CIF2, data, meanb, numBio = numBio, pREvec)
          Predraw2[j, jj] <- P2us
        }
      }
      
    } 
    
    
    for (jj in 1:N.ID) {
      Pred[[jj]] <- data.frame(u, Predraw1[jj, ], Predraw2[jj, ])
      colnames(Pred[[jj]]) <- c("times", "CIF1", "CIF2")
    }
  } 
  else{
    
    # ~~~~~~~~~~~~~~
    # Single Failure
    # ~~~~~~~~~~~~~~
    
    # intialize
    betaList <- object$betaList
    sigma <- object$sigma
    gamma1 <- object$gamma1
    
    alpha1flat <- object$alpha1
    
    if (latAsso == "sre") {
      
      alpha1 <- vector("list", numBio)
      index <- 1
      
      for (g in 1:numBio) {
        alpha1[[g]] <- alpha1flat[index:(index + pRE[g] - 1)]
        index <- index + pRE[g]
      }
      
    } else {
      
      # present/presentlp: one association parameter per biomarker
      alpha1 <- alpha1flat
      
    }
    
    H01 <- object$H01
    Sig <- object$Sig
    
    Predraw1 <- matrix(0, nrow = nrow(cnewdata2), ncol = length(u))
    lengthu <- length(u)
    
    subNDy <- YList <- XList <- y.obs <-  ZList <- vector("list", numBio)
    
    for (j in 1:N.ID) {
      
      stime <- as.numeric(Last.time[j])
      
      for(g in 1:numBio){
        
        # retain longitudinal at/before time s
        subNDy[[g]] <- ynewdata2[[g]][ynewdata2[[g]][, ID] == commonyID[j] &
                                        ynewdata2[[g]][, obs.time] <= stime, , drop = FALSE]
        
        y.obs[[g]][[j]] <- data.frame(ynewdata[ynewdata[, ID] == commonyID[j] &
                                                 ynewdata[, obs.time] <= stime, c(obs.time, Yvar[[g]][1])])
        
        # order longitudinal observations by time
        subNDy[[g]] <- subNDy[[g]][order(subNDy[[g]][, obs.time]), ,
                                   drop = FALSE]
        
        if (nrow(subNDy[[g]]) == 0) {
          stop("No longitudinal observations available by Last.time.")
        }
        
        YList[[g]] <- subNDy[[g]][, Yvar[[g]][1]]
        XList[[g]] <- as.matrix(data.frame(1, subNDy[[g]][, Yvar[[g]][-1]]))
        if (pRE[g] == 1) {
          ZList[[g]] <- matrix(1, ncol = 1, nrow = length(YList[[g]]))
        } else {
          ZList[[g]] <- as.matrix(data.frame(1, subNDy[[g]][, bvar1[[g]], drop = FALSE]))
        }
        
      }
      
      subNDc <- cnewdata2[cnewdata2[, ID] == commonyID[j], ]
      
      Xs_i <- NULL
      Zs_i <- NULL
      
      if (landmark && latAsso %in% c("present", "presentlp")) {
        
        Xs_i <- vector("list", numBio)
        Zs_i <- vector("list", numBio)
        
        for (g in 1:numBio) {
          
          # fixed-effects design for present value
          if (latAsso == "present") {
            
            Xs_i[[g]] <- XList[[g]][nrow(XList[[g]]), , drop = FALSE]
            
            time_fe <- which(Yvar[[g]][-1] == obs.time) + 1 # +1 for intercept
            
            if (length(time_fe) > 0) {
              # Xs_i[[g]][1, time_fe] <- lm.time
              Xs_i[[g]][1, time_fe] <- stime
            }
          }
          
          # random-effects design
          if (pRE[g] == 1) {
            
            Zs_i[[g]] <- matrix(1, nrow = 1, ncol = 1)
            
          } else {
            
            Zs_i[[g]] <- ZList[[g]][nrow(ZList[[g]]), , drop = FALSE]
            
            time_re <- which(bvar1[[g]] == obs.time) + 1
            
            if (length(time_re) > 0) {
              # Zs_i[[g]][1, time_re] <- lm.time
              Zs_i[[g]][1, time_re] <- stime
            }
          }
        }
      }
      
      # cumulative baseline hazard at prediction time
      CH01 <- CH(H01, stime)
      CH0u <- vector()
      for (jj in 1:lengthu) {
        CH0u[jj] <- CH(H01, u[jj])
      }
      
      W <- as.matrix(subNDc[1, Cvar[3:length(Cvar)]])
      
      data <- list(
        YList = YList,
        XList = XList,
        ZList = ZList,
        W = W,
        CH01 = CH01,
        beta = betaList,
        gamma1 = gamma1,
        alpha1 = alpha1,
        sigma = sigma,
        Sig = Sig,
        latAsso = latAsso,
        landmark = landmark,
        H01 = H01,
        stime = stime,
        Xs_i = Xs_i,
        Zs_i = Zs_i
      )
      
      opt <- optim(rep(0, nsig), logLikSFmv, data = data, method = "BFGS")
      meanb <- opt$par
      
      pREvec <- c()
      
      for(g in 1:numBio){
        pREvec[g] <- ncol(ZList[[g]])
      }
      
      # conditional surv prob at time u
      for (jj in 1:lengthu) {
        if (u[jj] > stime) {
          ## calculate the survival probabilities
          P1us <- P.us.mvjmcs(data, CH0u[jj], meanb, numBio, pREvec)
          Predraw1[j, jj] <- 1 - P1us
        } else {
          # set surv prob to 1 at or before time s
          Predraw1[j, jj] <- 1
        }
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
  
  
  if (min(u) < max(Last.time$Last.time)) {
    warning("Some prediction times precede Last.time; survival = 1, CIF = 0.")
  }
  
  sum
}


### PRINT WARNING ABOUT U < LAST.TIME



