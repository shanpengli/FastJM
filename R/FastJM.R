#' @useDynLib FastJM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom statmod  gauss.quad
#' @importFrom survival coxph Surv survfit
#' @importFrom dplyr left_join n across group_by row_number
#' @importFrom nlme lme getVarCov lmeControl
#' @importFrom graphics abline axis lines mtext panel.smooth par segments title
#' @importFrom MASS mvrnorm
#' @importFrom timeROC timeROC
#' @importFrom magrittr %>%
#' @importFrom stats as.formula median na.omit rexp rnorm optim pnorm qqline qqnorm quantile vcov model.matrix model.frame runif pchisq complete.cases rbeta
NULL