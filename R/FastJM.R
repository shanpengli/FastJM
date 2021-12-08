#' @useDynLib FastJM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom statmod  gauss.quad
#' @importFrom survival coxph
#' @importFrom dplyr left_join
#' @importFrom nlme lme getVarCov lmeControl
#' @importFrom mvtnorm rmvt
NULL