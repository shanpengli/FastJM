GetE <- function(beta, gamma1, gamma2, alpha1, alpha2, H01, H02, 
                 Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix,
                 method, Posbi, Poscov) {
  
  
  n <- nrow(X2)
  p1a <- ncol(Z)
  
  CUH01 <- rep(0, n)
  CUH02 <- rep(0, n)
  HAZ01 <- rep(0, n)
  HAZ02 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  CumuH02 <- cumsum(H02[, 3])
  
  getHazard(CumuH01, CumuH02, survtime, cmprsk, H01, H02, CUH01, CUH02, HAZ01, HAZ02)
  
  if (method == "standard") {
    status = getECstandard(beta, gamma1, gamma2, alpha1, alpha2,
                           Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, 
                   wsmatrix, CUH01, CUH02, HAZ01, HAZ02)
  } else {
    status = getECpseudo(beta, gamma1, gamma2, alpha1, alpha2,
                           Sig, sigma, Z, X1, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix,
                           wsmatrix, CUH01, CUH02, HAZ01, HAZ02, Posbi, Poscov)

  }
  
  
  return(status)
  
}

GetE.JMH <- function(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, 
                 Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix) {
  
  
  n <- nrow(X2)
  p1a <- ncol(Z)
  
  CUH01 <- rep(0, n)
  CUH02 <- rep(0, n)
  HAZ01 <- rep(0, n)
  HAZ02 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  CumuH02 <- cumsum(H02[, 3])
  
  getHazard(CumuH01, CumuH02, survtime, cmprsk, H01, H02, CUH01, CUH02, HAZ01, HAZ02)
  
  
  status = getEC(beta, tau, gamma1,  gamma2, alpha1, alpha2, vee1, vee2,  H01, 
                 H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, 
                 wsmatrix, CUH01, CUH02, HAZ01, HAZ02)
  
  return(status)
  
}

GetEad.JMH <- function(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, 
                       Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, 
                       wsmatrix, initial.optimizer, n_threads) {
  
  
  n <- nrow(X2)
  p1a <- ncol(Z)
  nsig <- p1a + 1
  CUH01 <- rep(0, n)
  CUH02 <- rep(0, n)
  HAZ01 <- rep(0, n)
  HAZ02 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  CumuH02 <- cumsum(H02[, 3])
  
  getHazard(CumuH01, CumuH02, survtime, cmprsk, H01, H02, CUH01, CUH02, HAZ01, HAZ02)
  
  subject_data <- lapply(seq_len(n), function(i) {
    rows <- if (i < n) {
      mdataS[i]:(mdataS[i + 1L] - 1L)
    } else {
      mdataS[i]:length(Y)
    }
    
    list(
      subject = i,
      Y        = Y[rows],
      X        = X1[rows, , drop = FALSE],
      Z        = Z[rows, , drop = FALSE],
      W        = W[rows, , drop = FALSE],
      X2       = X2[i, , drop = FALSE],
      CH01     = CUH01[i],
      CH02     = CUH02[i],
      HAZ01    = HAZ01[i],
      HAZ02    = HAZ02[i],
      D        = cmprsk[i]
    )
  })
  n_threads.opt <- max(floor(n_threads/2), 1)
  future::plan(
    future::multisession,
    workers = n_threads.opt
  )
  
  results <- future.apply::future_lapply(
    subject_data,
    function(data_i) {
      # Add smaller shared model parameters.
      i <- data_i$subject
      data_i$beta   <- beta
      data_i$tau    <- tau
      data_i$gamma1 <- gamma1
      data_i$gamma2 <- gamma2
      data_i$alpha1 <- alpha1
      data_i$alpha2 <- alpha2
      data_i$nu1    <- vee1
      data_i$nu2    <- vee2
      data_i$Sig    <- Sig
      
      opt <- optim(
        par = rep(0, nsig),
        fn = logLikCR.learn.JMH,
        data = data_i,
        method = initial.optimizer,
        hessian = TRUE
      )
      
      posterior_var <- tryCatch(
        solve(opt$hessian),
        error = function(e) {
          warning(sprintf(
            "Could not invert Hessian for subject %d: %s",
            i,
            conditionMessage(e)
          ))
          
          matrix(NA_real_, nsig, nsig)
        }
      )
      
      list(
        mode = opt$par,
        variance = posterior_var,
        convergence = opt$convergence,
        subject = i
      )
    },
    future.scheduling = 2
  )
  future::plan(future::sequential)
  posterior.mode <- do.call(
    rbind,
    lapply(results, `[[`, "mode")
  )
  posterior.var <- do.call(
    rbind,
    lapply(results, `[[`, "variance")
  )
  
  status = getECad(beta, tau, gamma1,  gamma2, alpha1, alpha2, vee1, vee2,  H01,
                   H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix,
                   wsmatrix, CUH01, CUH02, HAZ01, HAZ02, posterior.mode, posterior.var, n_threads)
  
  return(status)
  
}