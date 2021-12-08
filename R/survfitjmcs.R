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
##' @param M the number of Monte Carlo samples to be generated. Default is 200.
##' @param ... further arguments passed to or from other methods.
##' @export
##' 
survfitjmcs <- function(object, seed = 100, ynewdata = NULL, cnewdata = NULL, 
                           u = NULL, M = 200, ...) {
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
  if (!object$CompetingRisk) {
    H01 <- object$H01
    if (max(u) > H01[nrow(H01), 1])
      stop(paste("The current joint model cannot predict the conditional 
         survival probabilities later than the last observed time of the object. 
               The last observed time is", max(H01[, 1])))
  } else {
    H01 <- object$H01
    H02 <- object$H02
    if (max(u) > H01[nrow(H01), 1] | max(u) > H02[nrow(H02), 1])
      stop(paste("The current joint model cannot predict the conditional 
         survival probabilities later than the last observed time of the object. 
               The last observed time for risk 1 and 2 is", max(H01[, 1]), "and", max(H02[, 1])))
  }

  
  long.formula <- object$LongitudinalSubmodel
  Yvar <- all.vars(long.formula)
  surv.formula <- object$SurvivalSubmodel
  Cvar <- all.vars(surv.formula)
  bvar <- all.vars(object$random)
  if (!(bvar[length(bvar)] %in% colnames(ynewdata)))
    stop(paste("The ID variable", bvar[length(bvar)], "is not found in ynewdata."))
  if (!(bvar[length(bvar)] %in% colnames(cnewdata)))
    stop(paste("The ID variable", bvar[length(bvar)], "is not found in cnewdata."))
  
  ## dynamic prediction 
  ## Monte Carlo simulation
  # ND <- ydata2[ydata2$ID %in% c(23, 878, 230, 113), ]
  ID <- unique(ynewdata[, bvar[length(bvar)]])
  N.ID <- length(ID)
  cID <- cnewdata[, bvar[length(bvar)]]
  if (prod(ID == cID) == 0) {
    stop("The order of subjects in ydata doesn't match with cnewdata.")
  }
  set.seed(seed)
  nbeta <- length(object$beta)
  ngamma <- length(object$gamma1)
  nnu <- length(object$nu1)
  nsig <- nnu
  lengthu <- length(u)
  if (!object$CompetingRisk) {
    Psi <- c(object$beta, object$gamma1, object$nu1, object$sigma)
    for (l in 1:nsig) Psi <- c(Psi, object$Sig[l, l])
    if (nsig == 2) Psi <- c(Psi, object$Sig[1, 2])
    if (nsig == 3) {
      Psi <- c(Psi, object$Sig[1, 2])
      Psi <- c(Psi, object$Sig[2, 3])
      Psi <- c(Psi, object$Sig[1, 3])
    }
    covPsi <- vcov(object)
    Psi.MC <- mvrnorm(n = M, Psi, covPsi, tol = 1e-6, empirical = FALSE)
    Pred <- list()
    y.obs <- list()
    if (length(bvar) > 1) bvar1 <- bvar[1:(length(bvar) - 1)]
    for (j in 1:N.ID) {
      subNDy <- ynewdata[ynewdata[, bvar[length(bvar)]] == ID[j], ]
      subNDc <- cnewdata[cnewdata[, bvar[length(bvar)]] == ID[j], ]
      y.obs[[j]] <- data.frame(subNDy[, c(bvar[1], Yvar[1])])
      allPi <- matrix(0, ncol = length(u), nrow = M)
      CH0 <- CH(H01, subNDc[, Cvar[1]])
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
      for (i in 1:M) {
        ##1. draw Psi
        psil <- Psi.MC[i, ]
        betal <- psil[1:nbeta]
        gammal <- psil[(nbeta+1):(nbeta+ngamma)]
        nul <- psil[(nbeta+ngamma+1):(nbeta+ngamma+nnu)]
        sigmal <- psil[nbeta+ngamma+nnu+1]
        Sigl <- matrix(0, ncol = nsig, nrow = nsig)
        for (l in 1:nsig) Sigl[l, l] <- psil[nbeta+ngamma+nnu+1+l]
        if (nsig == 2) Sigl[1, 2] <- Sigl[2, 1] <- psil[nbeta+ngamma+nnu+1+nsig+1]
        if (nsig == 3) {
          Sigl[1, 2] <- Sigl[2, 1] <- psil[nbeta+ngamma+nnu+1+nsig+1]
          Sigl[2, 3] <- Sigl[3, 2] <- psil[nbeta+ngamma+nnu+1+nsig+2]
          Sigl[1, 3] <- Sigl[3, 1] <- psil[nbeta+ngamma+nnu+1+nsig+3]
        }
        data <- list(Y, X, Z, X2, CH0, betal, gammal, nul, sigmal, Sigl)
        names(data) <- c("Y", "X", "Z", "X2", "CH0", "beta", "gamma", "nu", "sigma", "Sig")
        opt <- optim(rep(0, nsig), logLik, data = data, method = "BFGS", hessian = TRUE)
        meanb <- opt$par
        varb <- solve(opt$hessian)
        bl <- as.vector(mvtnorm::rmvt(1, delta = meanb, sigma = varb, df = 4))
        for (jj in 1:lengthu) {
          Pi <- P.us(data, CH0u[jj], bl)
          allPi[i, jj] <- Pi
        }
      }
      allPi <- as.data.frame(allPi)
      colnames(allPi) <- u
      
      subCP <- as.data.frame(matrix(0, nrow = length(u), ncol = 5))
      colnames(subCP) <- c("times", "Mean", "Median", "95%Lower", "95%Upper")
      for (b in 1:length(u)) {
        subCP[b, 1] <- u[b]
        subCP[b, 2] <- 1 - mean(allPi[, b])
        subCP[b, 3] <- median(1 - allPi[, b])
        subCP[b, 4] <- quantile(1 - allPi[, b], probs = 0.025)
        subCP[b, 5] <- quantile(1 - allPi[, b], probs = 0.975)
      }
      Pred[[j]] <- subCP 
    }
    names(Pred) <- ID
    sum <- list()
    sum$Pred <- Pred
    class(sum) <- "survfitjmcs"
    sum$Last.time <- cnewdata[, c(bvar[length(bvar)], Cvar[1])]
    sum$M <- M
    names(y.obs) <- ID
    sum$y.obs <- y.obs
    sum
    #Pred
    #cat("\nPrediction of Conditional Probabilities of Event\n\tbased on", M, "Monte Carlo samples\n\n")
    #print(Pred)
  }

}