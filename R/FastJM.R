#' @useDynLib FastJM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom statmod  gauss.quad
#' @importFrom survival coxph Surv survfit
#' @importFrom dplyr left_join n across group_by row_number
#' @importFrom nlme lme getVarCov lmeControl
#' @importFrom graphics abline axis lines mtext panel.smooth par segments title legend panel.smooth plot par box
#' @importFrom MASS mvrnorm
#' @importFrom pec ipcw
#' @importFrom magrittr %>%
#' @importFrom stats sd as.formula fitted median na.omit rexp residuals delete.response terms rnorm optim pnorm qqline qqnorm quantile vcov model.matrix model.frame runif pchisq complete.cases rbeta
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom utils combn head read.table modifyList
#' @importFrom caret groupKFold
#' @importFrom rlang .data
#' @importFrom tidycmprsk cuminc
#' @importFrom ggpubr ggarrange
NULL