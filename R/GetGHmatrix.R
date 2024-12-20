GetGHmatrix <- function(quadpoint = quadpoint, Sig = Sig) {
  
  gq_vals <- statmod::gauss.quad(n = quadpoint, kind = "hermite")
  xs <- gq_vals$nodes
  ws <- gq_vals$weights
  
  nSig <- nrow(Sig)
  xsmatrix <- wsmatrix <- matrix(0, nrow = nSig, ncol = quadpoint^nSig)
  xsmatrix[nSig, ] <- rep(xs, quadpoint^(nSig-1))
  Total <- NULL
  for (i in 1:(nSig-1)) {
    for (j in 1:quadpoint) {
      sub <- rep(xs[j], quadpoint^i)
      Total <- c(Total, sub)
    }
    xsmatrix[nSig-i, ] <- rep(Total, quadpoint^(nSig-1-i))
    Total <- NULL
  }
  xsmatrix <- t(xsmatrix)
  
  wsmatrix[nSig, ] <- rep(ws, quadpoint^(nSig-1))
  Total <- NULL
  for (i in 1:(nSig-1)) {
    for (j in 1:quadpoint) {
      sub <- rep(ws[j], quadpoint^i)
      Total <- c(Total, sub)
    }
    wsmatrix[nSig-i, ] <- rep(Total, quadpoint^(nSig-1-i))
    Total <- NULL
  }
  wsmatrix <- t(wsmatrix)
  
  return(list(xsmatrix = xsmatrix, wsmatrix = wsmatrix))
}