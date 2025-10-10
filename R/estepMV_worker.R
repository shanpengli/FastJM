estepMV_worker<- function(j, data,pREtotal) {
  subCUH01 <- data$CUH01[j]
  subCUH02 <- data$CUH02[j]
  subHAZ01 <- data$HAZ01[j]
  subHAZ02 <- data$HAZ02[j]
  subcmprsk <- data$cmprsk[j]
  subW <- t(data$W[j, , drop = FALSE])
  
  subdata <- list(
    beta   = data$beta,  
    gamma1 = data$gamma1,
    gamma2 = data$gamma2,
    alpha  = data$alpha,
    sigma  = data$sigma,
    Z      = data$Z[[j]],
    X1     = data$X1[[j]],
    Y      = data$Y[[j]],
    Sig    = data$Sig,
    CUH01  = subCUH01,
    CUH02  = subCUH02,
    HAZ01  = subHAZ01,
    HAZ02  = subHAZ02,
    cmprsk = subcmprsk,
    W      = subW
  )
  
  opt <- optim(
    par     = rep(0, pREtotal),
    fn      = getbSig,
    gr      = getbSig_grad,
    data    = subdata,
    method  = "BFGS",
    hessian = TRUE
  )

  list(mode = opt$par, ccov = chol(solve(opt$hessian)))
}