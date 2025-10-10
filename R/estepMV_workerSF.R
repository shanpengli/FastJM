estepMV_workerSF<- function(j, data,pREtotal) {
  subCUH01 <- data$CUH01[j]
  subHAZ01 <- data$HAZ01[j]
  subcmprsk <- data$cmprsk[j]
  subW <- t(data$W[j, , drop = FALSE])
  
  subdata <- list(
    beta   = data$beta,  
    gamma1 = data$gamma1,
    alpha  = data$alpha,
    sigma  = data$sigma,
    Z      = data$Z[[j]],
    X1     = data$X1[[j]],
    Y      = data$Y[[j]],
    Sig    = data$Sig,
    CUH01  = subCUH01,
    HAZ01  = subHAZ01,
    cmprsk = subcmprsk,
    W      = subW
  )
  
  opt <- optim(
    par     = rep(0, pREtotal),
    fn      = getbSigSF,
    gr      = getbSig_gradSF,
    data    = subdata,
    method  = "BFGS",
    hessian = TRUE
  )
  
  list(mode = opt$par, ccov = chol(solve(opt$hessian)))
}