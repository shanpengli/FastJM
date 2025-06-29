% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jmcs.R
\name{jmcs}
\alias{jmcs}
\title{Joint modeling of longitudinal continuous data and competing risks}
\usage{
jmcs(
  ydata,
  cdata,
  long.formula,
  random = NULL,
  surv.formula,
  REML = TRUE,
  quadpoint = NULL,
  maxiter = 10000,
  print.para = FALSE,
  initial.para = NULL,
  survinitial = TRUE,
  tol = 1e-04,
  method = "pseudo-adaptive",
  opt = "nlminb"
)
}
\arguments{
\item{ydata}{a longitudinal data frame in long format.}

\item{cdata}{a survival data frame with competing risks or single failure.
Each subject has one data entry.}

\item{long.formula}{a formula object with the response variable and fixed effects covariates
to be included in the longitudinal sub-model.}

\item{random}{a one-sided formula object describing the random effects part of the longitudinal sub-model.
For example, fitting a random intercept model takes the form \code{ ~ 1|ID}.
Alternatively. Fitting a random intercept and slope model takes the form \code{~ x1 + ... + xn|ID}.}

\item{surv.formula}{a formula object with the survival time, event indicator, and the covariates
to be included in the survival sub-model.}

\item{REML}{a logic object that indicates the use of REML estimator. Default is TRUE.}

\item{quadpoint}{the number of pseudo-adaptive Gauss-Hermite quadrature points.
to be chosen for numerical integration. Default is 6 which produces stable estimates in most dataframes.}

\item{maxiter}{the maximum number of iterations of the EM algorithm that the function will perform. Default is 10000.}

\item{print.para}{Print detailed information of each iteration. Default is FALSE, i.e., not to print the iteration details.}

\item{initial.para}{a list of initialized parameters for EM iteration. Default is NULL.}

\item{survinitial}{Fit a Cox model to obtain initial values of the parameter estimates. Default is TRUE.}

\item{tol}{Tolerance parameter. Default is 0.0001.}

\item{method}{Method for proceeding numerical integration in the E-step. Default is pseudo-adaptive.}

\item{opt}{Optimization method to fit a linear mixed effects model, either \code{nlminb} (default) or \code{optim}.}
}
\value{
Object of class \code{jmcs} with elements
\item{beta}{the vector of fixed effects for the linear mixed effects model.} 
\item{gamma1}{the vector of fixed effects for type 1 failure for the survival model.}
\item{gamma2}{the vector of fixed effects for type 2 failure for the survival model. 
Valid only if \code{CompetingRisk = TRUE}.}
\item{nu1}{the vector of association parameter(s) for type 1 failure.}
\item{nu2}{the vector of association parameter(s) for type 2 failure. Valid only if \code{CompetingRisk = TRUE}.}
\item{H01}{the matrix that collects baseline hazards evaluated at each uncensored event time for type 1 failure. 
The first column denotes uncensored event times, the second column the number of events, and the third columns 
the hazards obtained by Breslow estimator.}
\item{H02}{the matrix that collects baseline hazards evaluated at each uncensored event time for type 2 failure. 
The data structure is the same as \code{H01}. Valid only if \code{CompetingRisk = TRUE}.}
\item{Sig}{the variance-covariance matrix of the random effects.}
\item{sigma}{the variance of the measurement error for the linear mixed effects model.}
\item{iter}{the total number of iterations until convergence.}
\item{convergence}{convergence identifier: 1 corresponds to successful convergence, 
whereas 0 to a problem (i.e., when 0, usually more iterations are required).}
\item{vcov}{the variance-covariance matrix of all the fixed effects for both models.}
\item{sebeta}{the standard error of \code{beta}.}
\item{segamma1}{the standard error of \code{gamma1}.}
\item{segamma2}{the standard error of \code{gamma2}. 
Valid only if \code{CompetingRisk = TRUE}.}
\item{senu1}{the standard error of \code{nu1}.}
\item{senu2}{the standard error of \code{nu2}. Valid only if \code{CompetingRisk = TRUE}.}
\item{seSig}{the vector of standard errors of covariance of random effects.}
\item{sesigma}{the standard error of variance of measurement error for the linear mixed effects model.}
\item{loglike}{the log-likelihood value.}
\item{fitted}{a list with the fitted values:
  \describe{
  \item{resid}{the vector of estimated residuals for the linear mixed effects model.} 
  \item{fitted}{the vector of fitted values for the linear mixed effects model.}
  \item{fittedmar}{the vector of marginal fitted values for the linear mixed effects model.}
  \item{residmar}{the vector of estimated marginal residuals for the linear mixed effects model.}
  }
}
\item{fittedSurv}{the estimated survival rate evaluated at each uncensored event time.}
\item{FUNB}{the estimated random effects for each subject.}
\item{CompetingRisk}{logical value; TRUE if a competing event are accounted for.}
\item{quadpoint}{the number of Gauss Hermite quadrature points used for numerical integration.}
\item{ydata}{the input longitudinal dataset for fitting a joint model.
It has been re-ordered in accordance with descending observation times in \code{cdata}.}
\item{cdata}{the input survival dataset for fitting a joint model.
It has been re-ordered in accordance with descending observation times.}
\item{PropEventType}{a frequency table of number of events.}
\item{LongitudinalSubmodel}{the component of the \code{long.formula}.}
\item{SurvivalSubmodel}{the component of the \code{surv.formula}.}
\item{random}{the component of the \code{random}.}
\item{call}{the matched call.}
\item{Quad.method}{the quadrature rule used for integration. 
If pseudo-adaptive quadrature rule is used, then return \code{pseudo-adaptive}. 
Otherwise return \code{standard}.}
\item{id}{the grouping vector for the longitudinal outcome.}
}
\description{
Joint modeling of longitudinal continuous data and competing risks
}
\examples{

require(FastJM)
require(survival)
# Load a simulated longitudinal dataset
data(ydata)
# Load a simulated survival dataset with two competing events
data(cdata)
\donttest{
# Fit a joint model
fit <- jmcs(ydata = ydata, cdata = cdata, 
            long.formula = response ~ time + gender + x1 + race, 
            surv.formula = Surv(surv, failure_type) ~ x1 + gender + x2 + race, 
            random =  ~ time| ID)
fit
# Extract the parameter estimates of longitudinal sub-model fixed effects
fixef(fit, process = "Longitudinal")
# Extract the parameter estimates of survival sub-model fixed effects
fixef(fit, process = "Event")
# Obtain the random effects estimates for first 6 subjects 
head(ranef(fit))
# Obtain the variance-covariance matrix of all parameter estimates 
vcov(fit)
# Obtain the result summaries of the joint model fit
summary(fit, process = "Longitudinal")
summary(fit, process = "Event")
# Prediction of cumulative incidence for competing risks data
# Predict the conditional probabilities for two patients who are alive (censored)
ND <- ydata[ydata$ID \%in\% c(419, 218), ]
ID <- unique(ND$ID)
NDc <- cdata[cdata$ID  \%in\% ID, ]
survfit <- survfitjmcs(fit, 
                       ynewdata = ND, 
                       cnewdata = NDc, 
                       u = seq(3, 4.8, by = 0.2), 
                       method = "GH",
                       obs.time = "time")
survfit
PE <- PEjmcs(fit, seed = 100, landmark.time = 3, horizon.time = c(3.6, 4, 4.4), 
             obs.time = "time", method = "GH", 
             quadpoint = NULL, maxiter = 1000, n.cv = 3, 
             survinitial = TRUE)
Brier <- summary(PE, error = "Brier")
Brier

MAEQ <- MAEQjmcs(fit, seed = 100, landmark.time = 3, horizon.time = c(3.6, 4, 4.4), 
                 obs.time = "time", method = "GH", 
                 quadpoint = NULL, maxiter = 1000, n.cv = 3, 
                 survinitial = TRUE)
APE <- summary(MAEQ, digits = 3)
APE

## evaluate prediction accuracy of fitted joint model using cross-validated mean AUC
AUC <- AUCjmcs(fit, seed = 100, landmark.time = 3, horizon.time = c(3.6, 4, 4.4),
               obs.time = "time", method = "GH",
               quadpoint = NULL, maxiter = 1000, n.cv = 3, metric = "AUC")
summary(AUC, digits = 3)

## or using cross-validated mean C-index
Cindex <- AUCjmcs(fit, seed = 100, landmark.time = 3, horizon.time = c(3.6, 4, 4.4),
               obs.time = "time", method = "GH",
               quadpoint = NULL, maxiter = 1000, n.cv = 3, metric = "Cindex")
summary(Cindex, digits = 3)

}

}
\seealso{
\code{\link{ranef}, \link{fixef}, \link{fitted.jmcs}, 
\link{residuals.jmcs}, \link{survfitjmcs}, \link{plot.jmcs},
\link{vcov.jmcs}}
}
\author{
Shanpeng Li \email{lishanpeng0913@ucla.edu}
}
