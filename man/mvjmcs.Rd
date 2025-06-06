% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvjmcs.R
\name{mvjmcs}
\alias{mvjmcs}
\title{Joint modeling of multivariate longitudinal and competing risks data}
\usage{
mvjmcs(
  ydata,
  cdata,
  long.formula,
  random = NULL,
  surv.formula,
  maxiter = 10000,
  opt = "nlminb",
  tol = 0.005,
  print.para = TRUE,
  initial.para = NULL,
  Gamma1Cal,
  Gamma2Cal
)
}
\arguments{
\item{ydata}{A longitudinal data frame in long format.}

\item{cdata}{A survival data frame with competing risks or single failure. Each subject has one data entry.}

\item{long.formula}{A list of formula objects specifying fixed effects for each longitudinal outcome.}

\item{random}{A formula or list of formulas describing random effects structures (e.g., \code{~ 1|ID}).}

\item{surv.formula}{A formula for the survival sub-model, including survival time and event indicator.}

\item{maxiter}{Maximum number of EM iterations. Default is 10000.}

\item{opt}{Optimization method for mixed model. Default is \code{"nlminb"}.}

\item{tol}{Convergence tolerance for EM algorithm. Default is 0.0001.}

\item{print.para}{Logical; if \code{TRUE}, prints parameter values at each iteration.}

\item{initial.para}{Optional list of initialized parameters. Default is \code{NULL}.}
}
\value{
A list containing:
\item{output}{EM algorithm output from final iteration}
\item{re}{Estimated random effects for each subject}
\item{sigi}{Estimated random effects covariance matrices (posterior) for each subject}
\item{beta}{Estimated fixed effects for longitudinal models}
\item{sigmaout}{Biomarker error variance estimates}
\item{gamma1}{Fixed effects for cause 1}
\item{gamma2}{Fixed effects for cause 2 (if \code{CompetingRisk = TRUE})}
\item{alpha1}{Association parameters for cause 1}
\item{alpha2}{Association parameters for cause 2}
\item{SEest}{Estimated standard errors of all parameters}
\item{runtime}{Elapsed run time}
\item{iter}{Number of EM iterations run}
}
\description{
Joint modeling of multivariate longitudinal continuous data and competing risks
}
\details{
Function fits a joint model for multiple longitudinal outcomes and competing risks using a fast EM algorithm.
}
\examples{


  require(FastJM)
  require(survival

  # Fit joint model with two biomarkers
  fit <-mvjmcs(ydata, cdata, long.formula = list(Y1 ~ X11 + X12 + time, Y2 ~ X11 + X12 + time),
                random = list(~time| ID, ~1|ID),
                surv.formula =Surv(survtime, cmprsk) ~ X21 + X22, maxiter = 50, opt = "nlminb", tol = 0.001)
  fit
  
  # Obtain the variance-variance matrix of all parameter estimates
  
}
\seealso{
\code{\link{jmcs}}, \code{\link{survfitjmcs}}, \code{\link{AUCjmcs}}, \code{\link{MAEQjmcs}}, \code{\link{PEjmcs}}
}
