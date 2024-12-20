##' @title Prediction in Joint Models
##' @name survPAjmcs
##' @aliases survPAjmcs
##' @description This function computes R square.
##' @param object an object inheriting from class \code{jmcs}.
##' @param ynewdata a data frame that contains the longitudinal and covariate information for the subjects 
##' for which prediction of survival probabilities is required.
##' @param cnewdata a data frame that contains the survival and covariate information for the subjects 
##' for which prediction of survival probabilities is required.
##' @param u a numeric vector of times for which prediction survival probabilities are to be computed.
##' @param s a numeric saclar. This specifies the known time at which each of 
##' the subjects in cnewdata was known to be alive.
##' @param obs.time a character string of specifying a longitudinal time variable in ynewdata.
##' @param LOCF a logical value to indicate whether the last-observation-carried-forward approach applies to prediction. 
##' If \code{TRUE}, then \code{LOCFcovariate} and \code{clongdata} must be specified to indicate 
##' which time-dependent survival covariates are included for dynamic prediction. Default is FALSE.
##' @param LOCFcovariate a vector of string with time-dependent survival covariates if \code{LOCF = TRUE}. Default is NULL.
##' @param clongdata a long format data frame where time-dependent survival covariates are incorporated. Default is NULL.
##' @param quadpoint number of quadrature points used for estimating conditional probabilities.
##' @param ... further arguments passed to or from other methods. 
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}}
##' @export
##' 
survPAjmcs <- function(object, ynewdata = NULL, cnewdata = NULL, 
                        u = NULL, s = NULL, obs.time = NULL, 
                        LOCF = FALSE, LOCFcovariate = NULL, clongdata = NULL,
                        quadpoint = NULL, ...) {
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  if (is.null(ynewdata))
    stop("New longitudinal data for dynamic prediction is needed.")
  if (is.null(cnewdata))
    stop("New longitudinal data for dynamic prediction is needed.")
  if (is.null(u)) 
    stop("Please specify the future time for dynamic prediction.")   
  if (!is.vector(u)) 
    stop("u must be vector typed.")
  if (is.null(quadpoint)) {
    quadpoint <- object$quadpoint
  }
  if (is.null(obs.time)) {
    stop("Please specify a vector that represents the time variable from ydatanew.")
  } else {
    if (!obs.time %in% colnames(ynewdata)) {
      stop(paste0(obs.time, " is not found in ynewdata."))
    }
  }
  
  bvar <- all.vars(object$random)
  ID <- bvar[length(bvar)]
  
  if (!(ID %in% colnames(ynewdata)))
    stop(paste("The ID variable", ID, "is not found in ynewdata."))
  if (!(ID %in% colnames(cnewdata)))
    stop(paste("The ID variable", ID, "is not found in cnewdata."))
  
  if (LOCF) {
    if (is.null(clongdata)) {
      stop("Please provide a long format data frame that includes all longitudinal covariates for dynamic prediction.\n")
    } else {
      survival <- all.vars(object$SurvivalSubmodel)
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
  
  ynewdata <- ynewdata[, colnames(object$ydata)]
  ynewdataraw <- ynewdata
  cnewdata <- cnewdata[, colnames(object$cdata)]
  
  if (LOCF) {
    clongdata <- clongdata[clongdata[, obs.time] <= s, ]
    clongdataLOCF <- clongdata %>%
      dplyr::group_by(dplyr::across(ID)) %>%
      dplyr::filter(row_number() == n())
    cnewdata[, LOCFcovariate] <- clongdataLOCF[, LOCFcovariate]
  }
  
  ydata2 <- rbind(object$ydata, ynewdata)
  cdata2 <- rbind(object$cdata, cnewdata)
  
  
  getdum <- getdummy(long.formula = object$LongitudinalSubmodel,
                     surv.formula = object$SurvivalSubmodel,
                     random = object$random, ydata = ydata2, cdata = cdata2)
  
  ydata2 <- getdum$ydata
  cdata2 <- getdum$cdata
  
  Yvar <- colnames(ydata2)[-1]
  Cvar <- colnames(cdata2)[-1]
  bvar <- all.vars(object$random)
  
  ny <- nrow(ynewdata)
  nc <- nrow(cnewdata)
  Ny <- nrow(ydata2)
  Nc <- nrow(cdata2)
  
  ynewdata <- ydata2[c((Ny-ny+1):Ny), ]
  cnewdata <- cdata2[c((Nc-nc+1):Nc), ]
  
  
  nsig <- nrow(object$Sig)
  
  getGH <- GetGHmatrix(quadpoint = quadpoint, Sig = object$Sig)
  
  xsmatrix <- getGH$xsmatrix
  wsmatrix <- getGH$wsmatrix
  
  if (length(bvar) > 1) bvar1 <- bvar[1:(length(bvar) - 1)]
  yID <- unique(ynewdata[, ID])
  N.ID <- length(yID)
  cID <- cnewdata[, ID]
  if (prod(yID == cID) == 0) {
    stop("The order of subjects in ydata doesn't match with cnewdata.")
  }
  
  if (sum(cnewdata[, Cvar[1]] < s) > 0)
    stop(paste0("Not all subjects are at risk at time ", s, "."))
  
  CompetingRisk <- object$CompetingRisk
  
  if (CompetingRisk) {
    
    Pred <- list(list(), list())
    beta <- object$beta
    sigma <- object$sigma
    gamma1 <- object$gamma1
    gamma2 <- object$gamma2
    nu1 <- object$nu1
    nu2 <- object$nu2
    H01 <- object$H01
    H02 <- object$H02
    Sig <- object$Sig
    
    y.obs <- list()
    lengthu <- length(u)
    CH01 <- CH(H01, s)
    CH02 <- CH(H02, s)
    
    ### extract all uncensored event times
    cnewdata.time1 <- sort(unique(cnewdata[cnewdata[, Cvar[2]] == 1, Cvar[1]]))
    cnewdata.time2 <- sort(unique(cnewdata[cnewdata[, Cvar[2]] == 2, Cvar[1]]))
    
    CIF1all <- matrix(0, nrow = length(cnewdata.time1), ncol = nrow(cnewdata) + 1)
    CIF2all <- matrix(0, nrow = length(cnewdata.time2), ncol = nrow(cnewdata) + 1)
    
    CIF1all[, 1] <- cnewdata.time1
    CIF2all[, 1] <- cnewdata.time2
    
    timecif1 <- H01[, 1][H01[, 1] >= s]
    timecif2 <- H02[, 1][H02[, 1] >= s]
    
    P1usAll <- P2usAll <- NULL
    
    P1usAll <- cbind(P1usAll, timecif1)
    P2usAll <- cbind(P2usAll, timecif2)
    
    for (j in 1:N.ID) {
      subNDy <- ynewdata[ynewdata[, ID] == yID[j], ]
      subNDc <- cnewdata[cnewdata[, ID] == yID[j], ]
      y.obs[[j]] <- data.frame(ynewdataraw[ynewdataraw[, ID] == yID[j], c(obs.time, Yvar[1])])
      
      Y <- subNDy[, Yvar[1]]
      X <- data.frame(1, subNDy[, Yvar[2:length(Yvar)]])
      X <- as.matrix(X)
      if (nsig == 1) {
        Z <- matrix(1, ncol = 1, nrow = length(Y))
      } else {
        Z <- data.frame(1, subNDy[, bvar1])
        Z <- as.matrix(Z)
      }
      X2 <- as.matrix(subNDc[1, Cvar[3:length(Cvar)]])
      
      ## find out E(bi)
      data <- list(Y, X, Z, X2, CH01, CH02, beta, gamma1, gamma2, nu1, nu2, sigma, Sig)
      names(data) <- c("Y", "X", "Z", "X2", "CH01", "CH02", "beta",
                       "gamma1", "gamma2", "nu1", "nu2",  "sigma", "Sig")
      opt <- optim(rep(0, nsig), logLikCR, data = data, method = "BFGS", hessian = TRUE)
      meanb <- opt$par
      Poscov <- solve(opt$hessian)
      
      CIF <- getECIFall(beta, sigma, gamma1, gamma2, nu1,
                        nu2, Sig, Z, X, Y, as.vector(X2), H01, H02,
                        xsmatrix, wsmatrix, CH01, CH02, s, timecif1, timecif2, meanb, Poscov)
      
      P1us <- CIF$CIF1
      P2us <- CIF$CIF2
      
      P1usAll <- cbind(P1usAll, P1us)
      P2usAll <- cbind(P2usAll, P2us)
      
    }
    
    for (i in 1:length(cnewdata.time1)) {
      logic <- cnewdata.time1[i] >= P1usAll[, 1]
      index <- max(which(logic))
      CIF1all[i, 2:(nrow(cnewdata) + 1)] <- P1usAll[index, 2:(nrow(cnewdata) + 1)]
    }
    CIF1all <- ifelse(is.na(CIF1all), 0, CIF1all)
    
    for (i in 1:length(cnewdata.time2)) {
      logic <- cnewdata.time2[i] >= P2usAll[, 1]
      index <- max(which(logic))
      CIF2all[i, 2:(nrow(cnewdata) + 1)] <- P2usAll[index, 2:(nrow(cnewdata) + 1)]
    }
    CIF2all <- ifelse(is.na(CIF2all), 0, CIF2all)
    
    for (i in 1:lengthu) {
      PAmeasure1 <- pam.censor.cr(ftime = cnewdata[, Cvar[1]], fstatus = cnewdata[, Cvar[2]],
                                  tau = u[i], pred.cif = CIF1all[, -1], time.cif = CIF1all[, 1],
                                  event.type = 1)
      Pred[[1]][[i]] <- PAmeasure1
      PAmeasure2 <- pam.censor.cr(ftime = cnewdata[, Cvar[1]], fstatus = cnewdata[, Cvar[2]],
                                  tau = u[i], pred.cif = CIF2all[, -1], time.cif = CIF2all[, 1],
                                  event.type = 2)
      Pred[[2]][[i]] <- PAmeasure2
    }
    names(Pred[[1]]) <- names(Pred[[2]]) <- u
    names(Pred) <- c("Risk1", "Risk2")
    
  } else {
    
    Predraw <- matrix(0, nrow = nrow(cnewdata), ncol = length(u))
    beta <- object$beta
    sigma <- object$sigma
    gamma <- object$gamma1
    nu <- object$nu1
    H01 <- object$H01
    Sig <- object$Sig
    
    lengthu <- length(u)
    
    for (j in 1:N.ID) {
      subNDy <- ynewdata[ynewdata[, ID] == yID[j], ]
      subNDc <- cnewdata[cnewdata[, ID] == yID[j], ]
      y.obs[[j]] <- data.frame(ynewdata[ynewdata[, ID] == yID[j], c(obs.time, Yvar[1])])
      
      CH0 <- CH(H01, Last.time[j])
      
      CH0u <- vector()
      for (jj in 1:lengthu) {
        CH0u[jj] <- CH(H01, u[jj])
      }
      Y <- subNDy[, Yvar[1]]
      X <- data.frame(1, subNDy[, Yvar[2:length(Yvar)]])
      X <- as.matrix(X)
      if (nsig == 1) {
        Z <- matrix(1, ncol = 1, nrow = length(Y))
      } else {
        Z <- data.frame(1, subNDy[, bvar1])
        Z <- as.matrix(Z)
      }
      X2 <- as.matrix(subNDc[, Cvar[3:length(Cvar)]])
      
      ## find out E(bi)
      data <- list(Y, X, Z, X2, CH0, beta, gamma, nu, sigma, Sig)
      names(data) <- c("Y", "X", "Z", "X2", "CH0", "beta", "gamma", "nu", "sigma", "Sig")
      opt <- optim(rep(0, nsig), logLik, data = data, method = "BFGS", hessian = TRUE)
      meanb <- opt$par
      Poscov <- solve(opt$hessian)
      
      if (method == "Laplace") {
        for (jj in 1:lengthu) {
          Pi <- P.us(data, CH0u[jj], meanb)
          Predraw[j, jj] <- 1 - Pi
        }
        quadpoint = NULL
      } else {
        for (jj in 1:lengthu) {
          Predraw[j, jj] <- getES(beta, sigma, gamma, nu, Sig, Z, X, Y, 
                                  as.vector(X2), xsmatrix, wsmatrix, CH0, CH0u[jj],
                                  meanb, Poscov)
        }
      }
      
    }
    for (jj in 1:N.ID) {
      Pred[[jj]] <- data.frame(u, Predraw[jj, ])
      colnames(Pred[[jj]]) <- c("times", "PredSurv")
    }
  }
  
  names(y.obs) <- yID
  s <- data.frame(cID, s)
  sum <- list(Pred = Pred, s = s, y.obs = y.obs, quadpoint = quadpoint,
              CompetingRisk = CompetingRisk)
  class(sum) <- "survPAjmcs"
  sum
  
}