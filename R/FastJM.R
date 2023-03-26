#' @useDynLib FastJM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom statmod  gauss.quad
#' @importFrom survival coxph survfit
#' @importFrom dplyr left_join
#' @importFrom nlme lme getVarCov lmeControl
#' @importFrom caret groupKFold 
#' @importFrom stats as.formula median optim pnorm vcov model.matrix model.frame runif pchisq complete.cases
NULL