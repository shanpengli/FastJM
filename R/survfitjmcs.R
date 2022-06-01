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
##' @param simulate logical; if \code{TRUE}, a Monte Carlo approach is used to estimate conditional probabilities. 
##' Otherwise, Gauss-Hermite quadrature rule is used for numerical integration to estimate instead. 
##' Default is \code{TRUE}.
##' @param quadpoint number of quadrature points used for estimating conditional probabilities when \code{simulate = FALSE}. Default is 20.
##' @param ... further arguments passed to or from other methods. 
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{jmcs}}
##' @export
##' 
survfitjmcs <- function(object, seed = 100, ynewdata = NULL, cnewdata = NULL, 
                        u = NULL, M = 200, simulate = TRUE, quadpoint = NULL, ...) {
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
  
  bvar <- all.vars(object$random)
  if (!(bvar[length(bvar)] %in% colnames(ynewdata)))
    stop(paste("The ID variable", bvar[length(bvar)], "is not found in ynewdata."))
  if (!(bvar[length(bvar)] %in% colnames(cnewdata)))
    stop(paste("The ID variable", bvar[length(bvar)], "is not found in cnewdata."))
  
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
  
  ## dynamic prediction 
  ## Monte Carlo simulation
  ID <- unique(ynewdata[, bvar[length(bvar)]])
  N.ID <- length(ID)
  cID <- cnewdata[, bvar[length(bvar)]]
  if (prod(ID == cID) == 0) {
    stop("The order of subjects in ydata doesn't match with cnewdata.")
  }
  
  if (simulate) {
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
          b.old <- meanb
          ## simulate new random effects estimates using the Metropolis Hastings algorithm
          propose.bl <- as.vector(mvtnorm::rmvt(1, delta = meanb, sigma = varb, df = 4))
          dmvt.old <- mvtnorm::dmvt(b.old, meanb, varb, df = 4, TRUE)
          dmvt.propose <- mvtnorm::dmvt(propose.bl, meanb, varb, df = 4, TRUE)
          logpost.old <- -logLik(data, b.old)
          logpost.propose <- -logLik(data, propose.bl)
          ratio <- min(exp(logpost.propose + dmvt.old - logpost.old - dmvt.propose), 1)
          if (runif(1) <= ratio) {
            bl = propose.bl
          } else {
            bl = b.old
          }
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
      sum$CompetingRisk <- FALSE
      sum$simulate <- simulate
      sum
      
    } else {
      
      Psi <- c(object$beta, object$gamma1, object$gamma2, 
               object$nu1, object$nu2, object$sigma)
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
        allPi1 <- matrix(0, ncol = length(u), nrow = M)
        allPi2 <- matrix(0, ncol = length(u), nrow = M)
        s <-  as.numeric(subNDc[1, Cvar[1]])
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
        X2 <- as.matrix(subNDc[, Cvar[3:length(Cvar)]])
        for (i in 1:M) {
          ##1. draw Psi
          psil <- Psi.MC[i, ]
          betal <- psil[1:nbeta]
          gammal1 <- psil[(nbeta+1):(nbeta+ngamma)]
          gammal2 <- psil[(nbeta+ngamma+1):(nbeta+2*ngamma)]
          nul1 <- psil[(nbeta+2*ngamma+1):(nbeta+2*ngamma+nnu)]
          nul2 <- psil[(nbeta+2*ngamma+nnu+1):(nbeta+2*ngamma+2*nnu)]
          sigmal <- psil[nbeta+2*ngamma+2*nnu+1]
          Sigl <- matrix(0, ncol = nsig, nrow = nsig)
          for (l in 1:nsig) Sigl[l, l] <- psil[nbeta+2*ngamma+2*nnu+1+l]
          if (nsig == 2) Sigl[1, 2] <- Sigl[2, 1] <- psil[nbeta+2*ngamma+2*nnu+1+nsig+1]
          if (nsig == 3) {
            Sigl[1, 2] <- Sigl[2, 1] <- psil[nbeta+2*ngamma+2*nnu+1+nsig+1]
            Sigl[2, 3] <- Sigl[3, 2] <- psil[nbeta+2*ngamma+2*nnu+1+nsig+2]
            Sigl[1, 3] <- Sigl[3, 1] <- psil[nbeta+2*ngamma+2*nnu+1+nsig+3]
          }
          data <- list(Y, X, Z, X2, CH01, CH02, betal, gammal1, gammal2, nul1, nul2, sigmal, Sigl)
          names(data) <- c("Y", "X", "Z", "X2", "CH01", "CH02", "beta", 
                           "gamma1", "gamma2", "nu1", "nu2",  "sigma", "Sig")
          opt <- optim(rep(0, nsig), logLikCR, data = data, method = "BFGS", hessian = TRUE)
          meanb <- opt$par
          varb <- solve(opt$hessian)
          b.old <- meanb
          ## simulate new random effects estimates using the Metropolis Hastings algorithm
          propose.bl <- as.vector(mvtnorm::rmvt(1, delta = meanb, sigma = varb, df = 4))
          dmvt.old <- mvtnorm::dmvt(b.old, meanb, varb, df = 4, TRUE)
          dmvt.propose <- mvtnorm::dmvt(propose.bl, meanb, varb, df = 4, TRUE)
          logpost.old <- -logLikCR(data, b.old)
          logpost.propose <- -logLikCR(data, propose.bl)
          ratio <- min(exp(logpost.propose + dmvt.old - logpost.old - dmvt.propose), 1)
          if (runif(1) <= ratio) {
            bl = propose.bl
          } else {
            bl = b.old
          }
          for (jj in 1:lengthu) {
            ## calculate the CIF
            CIF1 <- CIF1.CR(data, H01, H02, s, u[jj], bl)
            P1us <- Pk.us(CIF1, data, bl)
            allPi1[i, jj] <- P1us
            
            CIF2 <- CIF2.CR(data, H01, H02, s, u[jj], bl)
            P2us <- Pk.us(CIF2, data, bl)
            allPi2[i, jj] <- P2us
          }
        }
        allPi1 <- as.data.frame(allPi1)
        colnames(allPi1) <- u
        
        allPi2 <- as.data.frame(allPi2)
        colnames(allPi2) <- u
        
        subCP1 <- as.data.frame(matrix(0, nrow = length(u), ncol = 5))
        colnames(subCP1) <- c("times", "Mean", "Median", "95%Lower", "95%Upper")
        for (b in 1:length(u)) {
          subCP1[b, 1] <- u[b]
          subCP1[b, 2] <- mean(allPi1[, b])
          subCP1[b, 3] <- median(allPi1[, b])
          subCP1[b, 4] <- quantile(allPi1[, b], probs = 0.025)
          subCP1[b, 5] <- quantile(allPi1[, b], probs = 0.975)
        }
        
        subCP2 <- as.data.frame(matrix(0, nrow = length(u), ncol = 5))
        colnames(subCP2) <- c("times", "Mean", "Median", "95%Lower", "95%Upper")
        for (b in 1:length(u)) {
          subCP2[b, 1] <- u[b]
          subCP2[b, 2] <- mean(allPi2[, b])
          subCP2[b, 3] <- median(allPi2[, b])
          subCP2[b, 4] <- quantile(allPi2[, b], probs = 0.025)
          subCP2[b, 5] <- quantile(allPi2[, b], probs = 0.975)
        }
        
        subCP <- list(subCP1, subCP2)
        names(subCP) <- c("Cumulative incidence probabilities for type 1 failure",
                          "Cumulative incidence probabilities for type 2 failure")
        
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
      sum$CompetingRisk <- TRUE
      sum$simulate <- simulate
      sum
    }
  } else {
    if (is.null(quadpoint)) {
      quadpoint = 20 
    }
    gq_vals <- statmod::gauss.quad(n = quadpoint, kind = "hermite")
    xs <- gq_vals$nodes
    ws <- gq_vals$weights
    p1a <- ncol(object$Sig)
    
    if (p1a == 3) {
      xsmatrix <- matrix(0, nrow = 3, ncol = quadpoint^3)
      wsmatrix <- xsmatrix
      xsmatrix[3, ] <- rep(xs, quadpoint^2)
      Total <- NULL
      for (i in 1:quadpoint) {
        sub <- rep(xs[i], quadpoint)
        Total <- c(Total, sub)
      }
      xsmatrix[2, ] <- rep(Total, quadpoint)
      Total <- NULL
      for (i in 1:quadpoint) {
        sub <- rep(xs[i], quadpoint^2)
        Total <- c(Total, sub)
      }
      xsmatrix[1, ] <- Total
      xsmatrix <- t(xsmatrix)
      
      wsmatrix[3, ] <- rep(ws, quadpoint^2)
      Total <- NULL
      for (i in 1:quadpoint) {
        sub <- rep(ws[i], quadpoint)
        Total <- c(Total, sub)
      }
      wsmatrix[2, ] <- rep(Total, quadpoint)
      Total <- NULL
      for (i in 1:quadpoint) {
        sub <- rep(ws[i], quadpoint^2)
        Total <- c(Total, sub)
      }
      wsmatrix[1, ] <- Total
      wsmatrix <- t(wsmatrix)
    } else if (p1a == 2) {
      xsmatrix <- matrix(0, nrow = 2, ncol = quadpoint^2)
      wsmatrix <- xsmatrix
      xsmatrix[2, ] <- rep(xs, quadpoint)
      Total <- NULL
      for (i in 1:quadpoint) {
        sub <- rep(xs[i], quadpoint)
        Total <- c(Total, sub)
      }
      xsmatrix[1, ] <- Total
      xsmatrix <- t(xsmatrix)
      
      wsmatrix[2, ] <- rep(ws, quadpoint)
      Total <- NULL
      for (i in 1:quadpoint) {
        sub <- rep(ws[i], quadpoint)
        Total <- c(Total, sub)
      }
      wsmatrix[1, ] <- Total
      wsmatrix <- t(wsmatrix)
    } else {
      xsmatrix <- matrix(0, nrow = 1, ncol = quadpoint)
      wsmatrix <- xsmatrix
      xsmatrix[1, ] <- xs
      xsmatrix <- t(xsmatrix)
      wsmatrix[1, ] <- ws
      wsmatrix <- t(wsmatrix)
    }
    if (length(bvar) > 1) bvar1 <- bvar[1:(length(bvar) - 1)]
    if (!object$CompetingRisk) {
      Predraw <- matrix(0, nrow = nrow(cnewdata), ncol = length(u))
      beta <- object$beta
      sigma <- object$sigma
      gamma <- object$gamma1
      nu <- object$nu1
      H01 <- object$H01
      Sig <- object$Sig
      
      y.obs <- list()
      lengthu <- length(u)
      Last.time <- cnewdata[, c(bvar[length(bvar)], Cvar[1])]
      CompetingRisk <- object$CompetingRisk
      for (j in 1:N.ID) {
        subNDy <- ynewdata[ynewdata[, bvar[length(bvar)]] == ID[j], ]
        subNDc <- cnewdata[cnewdata[, bvar[length(bvar)]] == ID[j], ]
        y.obs[[j]] <- data.frame(subNDy[, c(bvar[1], Yvar[1])])
      }
      names(y.obs) <- ID
      Y <- ynewdata[, Yvar[1]]
      X <- data.frame(1, ynewdata[, Yvar[2:length(Yvar)]])
      X <- as.matrix(X)
      if (p1a == 1) {
        Z <- matrix(1, ncol = 1, nrow = length(Y))
      } else {
        Z <- data.frame(1, ynewdata[, bvar1])
        Z <- as.matrix(Z)
      }
      X2 <- as.matrix(cnewdata[, Cvar[3:length(Cvar)]])
      
      survtime <- cnewdata[, Cvar[1]]
      
      mdata <- as.data.frame(table(ynewdata[, bvar[length(bvar)]]))
      colnames(mdata) <- c("ID", "ni")
      mdata[, "ID"] <- as.character(mdata[, "ID"])
      n <- nrow(mdata)
      mdata <- as.vector(mdata$ni)
      mdataS <- rep(0, n) 
      mdataS[1] <- 1
      mdataCum <- cumsum(mdata)
      mdata2 <- mdata - 1
      mdataS[2:n] <- mdataCum[2:n] - mdata2[2:n]
      
      for (jj in 1:lengthu) {
        Pus <- getPus(beta, gamma, nu, Sig, sigma, H01, Z, X, Y, X2,
                      mdata, mdataS, survtime, xsmatrix, wsmatrix, u[jj])
        Predraw[, jj] <- Pus
      }
      Pred <- list()
      for (jj in 1:N.ID) {
        Pred[[jj]] <- data.frame(u, Predraw[jj, ])
        colnames(Pred[[jj]]) <- c("times", "PredSurv")
      }
      names(Pred) <- ID
      sum <- list(Pred = Pred, Last.time = Last.time, y.obs = y.obs, simulate = simulate, quadpoint = quadpoint)
      class(sum) <- "survfitjmcs"
      sum
    } else {
      beta <- object$beta
      sigma <- object$sigma
      gamma1 <- object$gamma1
      gamma2 <- object$gamma2
      nu1 <- object$nu1
      nu2 <- object$nu2
      H01 <- object$H01
      H02 <- object$H02
      Sig <- object$Sig
      
      Predraw1 <- matrix(0, nrow = nrow(cnewdata), ncol = length(u))
      Predraw2 <- matrix(0, nrow = nrow(cnewdata), ncol = length(u))
      y.obs <- list()
      lengthu <- length(u)
      Last.time <- cnewdata[, c(bvar[length(bvar)], Cvar[1])]
      CompetingRisk <- object$CompetingRisk
      for (j in 1:N.ID) {
        subNDy <- ynewdata[ynewdata[, bvar[length(bvar)]] == ID[j], ]
        subNDc <- cnewdata[cnewdata[, bvar[length(bvar)]] == ID[j], ]
        y.obs[[j]] <- data.frame(subNDy[, c(bvar[1], Yvar[1])])
      }
      names(y.obs) <- ID
      Y <- ynewdata[, Yvar[1]]
      X <- data.frame(1, ynewdata[, Yvar[2:length(Yvar)]])
      X <- as.matrix(X)
      if (p1a == 1) {
        Z <- matrix(1, ncol = 1, nrow = length(Y))
      } else {
        Z <- data.frame(1, ynewdata[, bvar1])
        Z <- as.matrix(Z)
      }
      X2 <- as.matrix(cnewdata[, Cvar[3:length(Cvar)]])
      
      survtime <- cnewdata[, Cvar[1]]
      
      mdata <- as.data.frame(table(ynewdata[, bvar[length(bvar)]]))
      colnames(mdata) <- c("ID", "ni")
      mdata[, "ID"] <- as.character(mdata[, "ID"])
      n <- nrow(mdata)
      mdata <- as.vector(mdata$ni)
      mdataS <- rep(0, n) 
      mdataS[1] <- 1
      mdataCum <- cumsum(mdata)
      mdata2 <- mdata - 1
      mdataS[2:n] <- mdataCum[2:n] - mdata2[2:n]
      
      for (jj in 1:lengthu) {
        Pus <- getPusCR(beta, gamma1, gamma2, nu1, nu2, Sig, sigma, H01, 
                        H02, Z, X, Y, X2, mdata, mdataS, survtime, 
                        xsmatrix, wsmatrix, u[jj])
        Predraw1[, jj] <- Pus[, 1]
        Predraw2[, jj] <- Pus[, 2]
      }
      Pred <- list()
      for (jj in 1:N.ID) {
        Pred[[jj]] <- data.frame(u, Predraw1[jj, ], Predraw2[jj, ])
        colnames(Pred[[jj]]) <- c("times", "CIF1", "CIF2")
      }
      names(Pred) <- ID
      sum <- list(Pred = Pred, Last.time = Last.time, y.obs = y.obs, simulate = simulate, quadpoint = quadpoint)
      class(sum) <- "survfitjmcs"
      sum
      
    }
  }
  
  
}
  


