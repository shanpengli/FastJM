#' Simulated longitudinal data
#'
#' @description The \code{ydata} data frame has 3067 rows and 4 columns.
#'
#' @format This data frame contains the following columns:
#'
#'   \describe{
#'
#'   \item{\code{ID}}{patient identifier.}
#'
#'   \item{\code{response}}{response variable.}
#'
#'   \item{\code{time}}{visit time.}
#'
#'   \item{\code{x1}}{treatment indicator. \code{0} denotes the placebo group and \code{1} the treatment group.}
#'
#'   }
#' @usage data(ydata)
#'
"ydata"

#' Simulated competing risks data
#'
#' @description The \code{cdata} data frame has 1000 rows and 5 columns.
#'
#' @format This data frame contains the following columns:
#'
#'   \describe{
#'
#'   \item{\code{ID}}{patient identifier.}
#
#'   \item{\code{surv}}{event time.}
#'
#'   \item{\code{failure_type}}{event indicator. \code{0} denotes censoring, \code{1} risk 1, 
#'   and \code{2} risk 2.}
#'
#'   \item{\code{x1}}{continuous variable.}
#'
#'   \item{\code{x2}}{treatment indicator. \code{0} denotes the placebo group and \code{1} the treatment group.}
#'
#'   }
#' @usage data(cdata)
#'
"cdata"
