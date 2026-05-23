#' Dynamic predictions from fitted joint models
#'
#' @description
#' \code{survfitJM()} is the recommended user-facing function for dynamic
#' prediction from fitted FastJM joint models. It dispatches automatically
#' according to the class of the fitted model object.
#'
#' @param object A fitted joint model object returned by \code{jmcs()},
#'   \code{JMMLSM()}, or \code{mvjmcs()}.
#' @param ... Additional arguments passed to the model-specific prediction
#'   method.
#'
#' @return An object containing dynamic prediction results. The exact structure
#'   depends on the fitted model class.
#'
#' @details
#' This function is an S3 generic. Depending on the class of \code{object}, it
#' dispatches to the corresponding model-specific prediction routine:
#' \code{survfitJM.jmcs()}, \code{survfitJM.JMMLSM()}, or
#' \code{survfitJM.mvjmcs()}.
#'
#' The model-specific functions \code{survfitjmcs()},
#' \code{survfitJMMLSM()}, and \code{survfitmvjmcs()} are retained as
#' lower-level functions for backward compatibility.
#'
#' @seealso \code{\link{jmcs}}, \code{\link{JMMLSM}}, \code{\link{mvjmcs}},
#'   \code{\link{survfitjmcs}}, \code{\link{survfitJMMLSM}},
#'   \code{\link{survfitmvjmcs}}
#'
#' @export
survfitJM <- function(object, ...) {
  UseMethod("survfitJM")
}