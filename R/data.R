#' Simulated longitudinal data
#'
#' @description The \code{ydata} data frame has 3067 rows and 6 columns.
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
#'    \item{\code{gender}}{gender indicator.}
#'    
#'    \item{\code{race}}{race indicator.}
#'   }
#' @usage data(ydata)
#'
"ydata"

#' Simulated competing risks data correlated with ydata
#'
#' @description The \code{cdata} data frame has 1000 rows and 7 columns.
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
#'   \item{\code{gender}}{gender indicator.}
#'    
#'   \item{\code{race}}{race indicator.}
#'   }
#' @usage data(cdata)
#'
"cdata"

#' Simulated bivariate longitudinal data
#' @description The \code{mvydata} data frame has 4060 rows and 6 columns.
#'
#' @format This data frame contains the following columns:
#'
#'   \describe{
#'
#'   \item{\code{ID}}{patient identifier.}
#'   
#'   \item{\code{time}}{visit time.}
#'   
#'   \item{\code{Y1}}{response variable of biomarker 1.}
#'   
#'   \item{\code{Y2}}{response variable of biomarker 2.}
#'
#'   \item{\code{X11}}{X11.}
#'    
#'    \item{\code{X12}}{X12.}
#'    
#'   }
#' @usage data(mvydata)
#'
"mvydata"

#' Simulated competing risks data correlated with mvydata
#'
#' @description The \code{mvcdata} data frame has 500 rows and 5 columns.
#'
#' @format This data frame contains the following columns:
#'
#'   \describe{
#'
#'   \item{\code{ID}}{patient identifier.}
#
#'   \item{\code{survtime}}{event time.}
#'
#'   \item{\code{cmprsk}}{event indicator. \code{0} denotes censoring, \code{1} risk 1, 
#'   and \code{2} risk 2.}
#'   
#'   \item{\code{X21}}{X21.}
#'    
#'   \item{\code{X22}}{X22.}
#'   }
#' @usage data(mvcdata)
#'
"mvcdata"

#' Simulated longitudinal data with within-subject variance
#'
#' @description The \code{ydatah} data frame has 1353 rows and 6 columns.
#'
#' @format This data frame contains the following columns:
#'
#'   \describe{
#'
#'   \item{\code{ID}}{patient identifier.}
#'
#'   \item{\code{Y}}{response variable.}
#'
#'   \item{\code{time}}{visit time.}
#'
#'   \item{\code{Z1}}{treatment indicator. \code{0} denotes the placebo group and \code{1} the treatment group.}
#'    
#'    \item{\code{Z2}}{continuous variable..}
#'    
#'    \item{\code{Z3}}{continuous variable..}
#'   }
#' @usage data(ydatah)
#'
"ydatah"

#' Simulated competing risks data where event hazards depend on within-subject variance
#'
#' @description The \code{cdatah} data frame has 200 rows and 6 columns.
#'
#' @format This data frame contains the following columns:
#'
#'   \describe{
#'
#'   \item{\code{ID}}{patient identifier.}
#
#'   \item{\code{survtime}}{event time.}
#'
#'   \item{\code{cmprsk}}{event indicator. \code{0} denotes censoring, \code{1} risk 1, 
#'   and \code{2} risk 2.}
#'
#'   \item{\code{var1}}{treatment indicator. \code{0} denotes the placebo group and \code{1} the treatment group.}
#'
#'   \item{\code{var2}}{continuous variable.}
#'
#'   \item{\code{var3}}{continuous variable.}
#'   }
#' @usage data(cdatah)
#'
"cdatah"

