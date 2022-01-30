#' @useDynLib FastJM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom statmod  gauss.quad
#' @importFrom survival coxph
#' @importFrom dplyr left_join
#' @importFrom nlme lme getVarCov lmeControl
#' @importFrom mvtnorm rmvt
#' @importFrom graphics abline axis lines mtext panel.smooth par segments
#' @importFrom MASS mvrnorm
#' @importFrom stats as.formula median optim pnorm qqline qqnorm quantile vcov model.matrix
NULL