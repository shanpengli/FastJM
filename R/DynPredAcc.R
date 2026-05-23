##' Dynamic prediction accuracy metrics for fitted joint models
##'
##' Computes dynamic prediction accuracy measures for fitted joint models of
##' class \code{jmcs}, \code{JMMLSM}, or \code{mvjmcs} using grouped
##' cross-validation.
##'
##' At each cross-validation fold, the function refits the joint model on the
##' training set, obtains subject-specific dynamic predictions for validation
##' subjects at the requested horizon times, and evaluates the requested
##' prediction accuracy metrics.
##'
##' Supported metrics include:
##' \describe{
##'   \item{\code{"AUC"}}{Time-dependent area under the ROC curve.}
##'   \item{\code{"Cindex"}}{Concordance index for dynamic predictions.}
##'   \item{\code{"Brier Score"}}{Inverse-probability-of-censoring weighted
##'   Brier score.}
##'   \item{\code{"MAE"}}{Inverse-probability-of-censoring weighted mean
##'   absolute error.}
##'   \item{\code{"MAEQ"}}{Quantile-based calibration summaries comparing
##'   empirical and predicted cumulative incidence rates or survival probabilities.}
##' }
##'
##' For competing-risks models, event-specific cumulative incidence predictions
##' are evaluated. For single-event models, conditional survival probabilities
##' are evaluated.
##'
##' @title Dynamic prediction accuracy metrics for joint models
##' @name DynPredAcc
##' @aliases DynPredAcc
##'
##' @param seed A numeric value used to set the random seed for cross-validation.
##' Default is \code{100}.
##' @param object A fitted object of class \code{jmcs}, \code{JMMLSM}, or
##' \code{mvjmcs}.
##' @param landmark.time A numeric value specifying the landmark time at which
##' dynamic prediction begins.
##' @param horizon.time A numeric vector of future times at which predicted
##' probabilities are evaluated.
##' @param obs.time A character string specifying the longitudinal time variable
##' in the longitudinal data.
##' @param method A character string specifying the approximation method used for
##' dynamic prediction. Available options are \code{"Laplace"} and \code{"GH"}.
##' This argument is not used for objects of class \code{mvjmcs}.
##' @param quadpoint Number of Gauss--Hermite quadrature points used when
##' \code{method = "GH"}. If \code{NULL}, the value stored in \code{object} is
##' used.
##' @param maxiter Maximum number of EM iterations allowed when refitting the
##' model within each cross-validation fold. If \code{NULL}, the model-specific
##' default is used.
##' @param n.cv Number of cross-validation folds. Default is \code{3}.
##' @param quantile.width Numeric value specifying the width of the quantile
##' groups used for \code{"MAEQ"} summaries. Default is \code{0.25}. The
##' reciprocal of \code{quantile.width} must be an integer.
##' @param opt Optimization method used when refitting the model in each fold.
##' Available options are \code{"nlminb"} and \code{"optim"}.
##' @param LOCF Logical value indicating whether the last-observation-carried-
##' forward approach is applied for time-dependent survival covariates during
##' prediction. Default is \code{FALSE}.
##' @param LOCFcovariate A character vector specifying the time-dependent
##' survival covariates to be updated by last observation carried forward when
##' \code{LOCF = TRUE}. Default is \code{NULL}.
##' @param clongdata A long-format data frame containing time-dependent survival
##' covariates used when \code{LOCF = TRUE}. Default is \code{NULL}.
##' @param metrics A character vector specifying which prediction accuracy
##' metrics to compute. Available options are \code{"AUC"}, \code{"Cindex"},
##' \code{"Brier Score"}, \code{"MAE"}, and \code{"MAEQ"}. Default is
##' \code{c("AUC", "Cindex", "Brier Score", "MAE", "MAEQ")}.
##' @param ... Further arguments passed to model-specific methods.
##'
##' @return
##' An object of class \code{DynPredAcc}, returned as a list containing:
##' \describe{
##'   \item{\code{jm.class}}{The class of the fitted joint model.}
##'   \item{\code{n.cv}}{The number of cross-validation folds.}
##'   \item{\code{landmark.time}}{The landmark time used for dynamic prediction.}
##'   \item{\code{horizon.time}}{The vector of horizon times used for evaluation.}
##'   \item{\code{method}}{The dynamic prediction approximation method used, if applicable.}
##'   \item{\code{quadpoint}}{The number of quadrature points used, if applicable.}
##'   \item{\code{CompetingRisk}}{Logical value indicating whether the fitted model
##'   accounts for competing risks.}
##'   \item{\code{seed}}{The random seed used for cross-validation.}
##'   \item{\code{metrics}}{The requested evaluation metrics.}
##'   \item{\code{quantile.width}}{The width of quantile groups used for
##'   quantile-based calibration summaries.}
##'   \item{\code{AUC.cv}}{A list of fold-specific time-dependent AUC estimates,
##'   returned when \code{"AUC"} is requested.}
##'   \item{\code{Cindex.cv}}{A list of fold-specific concordance index estimates,
##'   returned when \code{"Cindex"} is requested.}
##'   \item{\code{Brier.cv}}{A list of fold-specific inverse-probability-of-
##'   censoring weighted Brier scores, returned when \code{"Brier Score"} is
##'   requested.}
##'   \item{\code{MAE.cv}}{A list of fold-specific inverse-probability-of-
##'   censoring weighted mean absolute errors, returned when \code{"MAE"} is
##'   requested.}
##'   \item{\code{MAEQ.cv}}{A list of fold-specific quantile-based calibration
##'   summaries, returned when \code{"MAEQ"} is requested.}
##' }
##'
##' @details
##' The function performs grouped cross-validation at the subject level. Within
##' each fold, the model is refit on the training data, and dynamic predictions
##' are computed for validation subjects who remain under observation beyond the
##' landmark time and who have longitudinal observations observed up to that
##' time. For single-event models, prediction accuracy is based on conditional
##' survival probabilities. For competing-risks models, prediction accuracy is
##' based on event-specific cumulative incidence functions.
##'
##' @seealso \code{\link{jmcs}}, \code{\link{JMMLSM}}, \code{\link{mvjmcs}},
##' \code{\link{survfitjmcs}}, \code{\link{survfitJMMLSM}}, \code{\link{survfitmvjmcs}}
##'
##' @export
##' 
DynPredAcc <- function(seed = 100, object, landmark.time = NULL, horizon.time = NULL,
                       obs.time = NULL, method = c("Laplace", "GH"),
                       quadpoint = NULL, maxiter = NULL, n.cv = 3,
                       quantile.width = 0.25,
                       opt = c("nlminb", "optim"),
                       LOCF = FALSE, LOCFcovariate = NULL, clongdata = NULL,
                       metrics = c("AUC", "Cindex", "Brier Score", "MAE", "MAEQ"), 
                       cpu.cores = 1, ...) {
  
  if (!any(class(object) %in% c("jmcs", "mvjmcs", "JMMLSM"))) {
    stop("Use only with 'jmcs', 'mvjmcs', or 'JMMLSM' objects.\n")
  }
  
  method <- match.arg(method)
  opt <- match.arg(opt)
  
  if (inherits(object, "jmcs")) {
    res <- DynPredAcc.jmcs(seed = seed, object = object, 
                           landmark.time = landmark.time, 
                           horizon.time = horizon.time,
                           obs.time = obs.time, method = method,
                           quadpoint = quadpoint, maxiter = maxiter, 
                           n.cv = n.cv,
                           opt = opt,
                           quantile.width = quantile.width,
                           LOCF = LOCF, LOCFcovariate = LOCFcovariate, 
                           clongdata = clongdata,
                           metrics = metrics)
  } else if (inherits(object, "JMMLSM")) {
    res <- DynPredAcc.JMMLSM(seed = seed, object = object, 
                             landmark.time = landmark.time, 
                             horizon.time = horizon.time,
                             obs.time = obs.time, method = method,
                             quadpoint = quadpoint, maxiter = maxiter, 
                             n.cv = n.cv,
                             opt = opt,
                             quantile.width = quantile.width,
                             LOCF = LOCF, LOCFcovariate = LOCFcovariate, 
                             clongdata = clongdata,
                             metrics = metrics)
  } else {
    res <- DynPredAcc.mvjmcs(seed = seed, object = object, 
                             landmark.time = landmark.time, 
                             horizon.time = horizon.time,
                             obs.time = obs.time,
                             maxiter = maxiter, 
                             n.cv = n.cv,
                             opt = opt,
                             quantile.width = quantile.width,
                             LOCF = LOCF,
                             LOCFcovariate = LOCFcovariate,
                             clongdata = clongdata,
                             metrics = metrics,
                             cpu.cores = cpu.cores)
    
  }
  return(res)
}