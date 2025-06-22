##' @title Prediction in Joint Models
##' @name survfitjmcs
##' @aliases survfitjmcs
##' @description This function computes the conditional probability of 
##' surviving later times than the last observed time for which a longitudinal 
##' measurement was available.
##' @param object an object inheriting from class \code{jmcs}.
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
##' @param method a character string specifying the type of probability approximation; if \code{Laplace}, then a first order estimator is computed.
##' If \code{GH}, then the standard Gauss-Hermite quadrature is used instead.
##' @param quadpoint number of quadrature points used for estimating conditional probabilities 
##' when \code{method = "GH"}. Default is NULL. If \code{method = "GH"}, 
##' then use the same amount of quadrature points obtained from \code{object}.
##' @param ... further arguments passed to or from other methods. 
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}}
##' @export
##' 
survfitjmcs <- function(object, seed = 100, ynewdata = NULL, cnewdata = NULL, 
                        u = NULL, Last.time = NULL, obs.time = NULL, 
                        LOCF = FALSE, LOCFcovariate = NULL, clongdata = NULL,
                        method = c("Laplace", "GH"), 
                        quadpoint = NULL, ...) {
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  if (is.null(ynewdata))
    stop("New longitudinal data for dynamic prediction is needed.")
  if (is.null(cnewdata))
    stop("New longitudinal data for dynamic prediction is needed.")
  if (is.null(u)) 
    stop("Please specify the future time for dynamic prediction.")   
  if (!method %in% c("Laplace", "GH"))
    stop("Please specify a method for probability approximation: Laplace or GH.")
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
  cnewdata <- cnewdata[, colnames(object$cdata)]
  
  if (LOCF) {
    clongdata <- clongdata[clongdata[, obs.time] <= Last.time, ]
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
  
  ynewdata2 <- ydata2[c((Ny-ny+1):Ny), ]
  cnewdata2 <- cdata2[c((Nc-nc+1):Nc), ]
  
  
  nsig <- nrow(object$Sig)
  
  getGH <- GetGHmatrix(quadpoint = quadpoint, Sig = object$Sig)
  
  xsmatrix <- getGH$xsmatrix
  wsmatrix <- getGH$wsmatrix
  
  if (length(bvar) > 1) bvar1 <- bvar[1:(length(bvar) - 1)]
  yID <- unique(ynewdata2[, ID])
  N.ID <- length(yID)
  cID <- cnewdata2[, ID]
  if (prod(yID == cID) == 0) {
    stop("The order of subjects in ydata doesn't match with cnewdata.")
  }
  
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
      } 
      if (length(Last.time) < nrow(cnewdata)) {
        stop("The last.time vector does not match cnewdata.")
      }
    }
  } else {
    Last.time <- cnewdata[, Cvar[1]]
  }
  
  Pred <- list()
  y.obs <- list()
  CompetingRisk <- object$CompetingRisk
  
  if (CompetingRisk) {

    beta <- object$beta
    sigma <- object$sigma
    gamma1 <- object$gamma1
    gamma2 <- object$gamma2
    nu1 <- object$nu1
    nu2 <- object$nu2
    H01 <- object$H01
    H02 <- object$H02
    Sig <- object$Sig

    Predraw1 <- matrix(0, nrow = nrow(cnewdata2), ncol = length(u))
    Predraw2 <- matrix(0, nrow = nrow(cnewdata2), ncol = length(u))
    lengthu <- length(u)

    for (j in 1:N.ID) {
      subNDy <- ynewdata2[ynewdata2[, ID] == yID[j], ]
      subNDc <- cnewdata2[cnewdata2[, ID] == yID[j], ]
      y.obs[[j]] <- data.frame(ynewdata[ynewdata[, ID] == yID[j], c(obs.time, Yvar[1])])

      s <-  as.numeric(Last.time[j])
      CH01 <- CH(H01, s)
      CH02 <- CH(H02, s)
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

      if (method == "GH") {
        for (jj in 1:lengthu) {
          ## calculate the CIF
          CIF <- getECIF(beta, sigma, gamma1, gamma2, nu1, nu2, Sig, Z, X, Y,
                         as.vector(X2), H01, H02,
                         xsmatrix, wsmatrix, CH01, CH02, s, u[jj], meanb, Poscov)
          
          P1us <- CIF$CIF1
          P2us <- CIF$CIF2

          Predraw1[j, jj] <- P1us
          Predraw2[j, jj] <- P2us
        }
      } else {

        for (jj in 1:lengthu) {
          ## calculate the CIF
          CIF1 <- CIF1.CR(data, H01, H02, s, u[jj], meanb)
          P1us <- Pk.us(CIF1, data, meanb)
          Predraw1[j, jj] <- P1us

          CIF2 <- CIF2.CR(data, H01, H02, s, u[jj], meanb)
          P2us <- Pk.us(CIF2, data, meanb)
          Predraw2[j, jj] <- P2us
        }
        quadpoint = NULL
      }

      for (jj in 1:N.ID) {
        Pred[[jj]] <- data.frame(u, Predraw1[jj, ], Predraw2[jj, ])
        colnames(Pred[[jj]]) <- c("times", "CIF1", "CIF2")
      }

    }

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
      subNDy <- ynewdata2[ynewdata2[, ID] == yID[j], ]
      subNDc <- cnewdata2[cnewdata2[, ID] == yID[j], ]
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
  
  names(y.obs) <- names(Pred) <- yID
  Last.time <- data.frame(cID, Last.time)
  sum <- list(Pred = Pred, Last.time = Last.time, y.obs = y.obs, method = method, quadpoint = quadpoint,
              CompetingRisk = CompetingRisk)
  class(sum) <- "survfitjmcs"
  sum
  
}

