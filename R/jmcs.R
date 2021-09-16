##' Joint modeling of longitudinal continuous data and competing risks
##' @title Joint Modelling for Continuous outcomes
##' @param ydata a longitudinal data frame in long format.
##' @param cdata a survival data frame with competing risks or single failure.
##' Each subject has one data entry.
##' @param long.formula a formula object with the response variable and fixed effects covariates
##' to be included in the longitudinal sub-model.
##' @param random a one-sided formula object describing the random effects part of the longitudinal sub-model.
##' For example, fitting a random intercept model takes the form ~ 1|ID.
##' Alternatively. Fitting a random intercept and slope model takes the form ~ x1 + ... + xn|ID.
##' @param surv.formula a formula object with the survival time, event indicator, and the covariates
##' to be included in the survival sub-model.
##' @param REML a logic object that indicates the use of REML estimator. Default is TRUE.
##' @param point the number of pseudo-adaptive Gauss-Hermite quadrature points
##' to be chosen for numerical integration. Default is 6 which produces stable estimates in most dataframes.
##' @param maxiter the maximum number of iterations of the EM algorithm that the function will perform. Default is 10000.
##' @param do.trace Print detailed information of each iteration. Default is FALSE, i.e., not to print the iteration details.
##' @param survinitial Fit a Cox model to obtain initial values of the parameter estimates. Default is TRUE.
##' @param tol Tolerance parameter. Default is 0.001.
##' @return Object of class \code{FastJM} with elements
##'   \tabular{ll}{
##'       \code{vcmatrix}    \tab  The variance-covariance matrix for all the parameters. The parameters are in the order: \eqn{\beta}, \eqn{\sigma^2}, \eqn{\gamma}, \eqn{\nu}, and \eqn{\Sigma}. The elements in \eqn{\Sigma} are output in the order along the main diagonal line, then the second main diagonal line, and so on. \cr
##'       \code{betas} \tab The point  estimates of \eqn{\beta}. \cr
##'       \code{se_betas} \tab The standard error estimate of \eqn{\beta}. \cr
##'       \code{gamma_matrix} \tab  The point  estimate of \eqn{\gamma}. \cr
##'       \code{se_gamma_matrix}   \tab  The standard error estimate of \eqn{\gamma}. \cr
##'       \code{vee1_estimate} \tab The point  estimate of \eqn{\nu_1}. \cr
##'       \code{se_vee1_estimate}    \tab The standard error estimate of \eqn{\nu_1}. \cr
##'       \code{vee2_estimate} \tab The point  estimate of \eqn{\nu_2}. \cr
##'       \code{se_vee2_estimate}    \tab The standard error estimate of \eqn{\nu_2}. \cr
##'       \code{sigma2_val}     \tab  The point estimate of \eqn{\sigma^2}.\cr
##'       \code{se_sigma2_val}     \tab  The standard error estimate of \eqn{\sigma^2}.\cr
##'       \code{sigma_matrix}     \tab The point estimate of \eqn{\Sigma} (only the upper triangle portion of the matrix is output).\cr
##'       \code{se_sigma}     \tab The standard error estimate of \eqn{\Sigma}.The standard errors are given in this order: main diagonal, the second main diagonal, and so on. \cr
##'       \code{loglike}     \tab Log Likelihood.\cr
##'   }
##' @export
##'

jmcs <- function(ydata, cdata, long.formula, random = NULL, surv.formula, REML = TRUE,
                 point = 6, maxiter = 10000, do.trace = FALSE, survinital = TRUE, tol = 0.0001)
{
  if (do.trace) {
    trace=1;
  }else{
    trace=0;
  }

  #Gaussian-Hermite quadrature nodes and weights
  #The dimension of xs/ws is half of the point value since they are symmetric

  if (point %% 2 == 1)
  {
    stop("Number of quadrature points can only be even!")
  }

  gq_vals <- statmod::gauss.quad(n = point, kind = "hermite")

  xs <- gq_vals$nodes[(point / 2 + 1) : point]

  ws <- gq_vals$weights[(point / 2 + 1) : point]

  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
  random.form <- all.vars(random)
  ID <- random.form[length(random.form)]
  if (length(random.form) == 1) {
    RE <- NULL
    model <- "intercept"
  } else {
    RE <- random.form[-length(random.form)]
    model <- "interslope"
  }

  cnames=colnames(cdata)
  ynames=colnames(ydata)

  ##variable check
  if (prod(long %in% ynames) == 0) {
    Fakename <- which(long %in% ynames == FALSE)
    stop(paste0("The variable ", long[Fakename], " not found"))
  }
  if (prod(survival %in% cnames) == 0) {
    Fakename <- which(survival %in% cnames == FALSE)
    stop(paste0("The variable ", survival[Fakename], " not found"))
  }
  if (!(ID %in% ynames)) {
    stop(paste0("ID column ", ID, " not found in long_data"))
  }
  if (!(ID %in% cnames)) {
    stop(paste0("ID column ", ID, " not found in surv_data"))
  }
  if (is.null(RE) & model == "interslope") {
    stop("Random effects covariates must be specified.")
  }
  yID <- unique(ydata[, ID])
  cID <- cdata[, ID]
  if (prod(yID == cID) == 0) {
    stop("The order of subjects in ydata doesn't match with cdata.")
  }


  ydim = dim(ydata)
  cdim = dim(cdata)
  mdata <- as.data.frame(table(ydata[, ID]))
  mdata <- as.data.frame(mdata[, 2])
  n1 = ydim[1]


  if (nrow(mdata) != cdim[1]) {
    stop(paste("The number of subjects in cdata are not consistent with ydata.
         Please make sure the cdata has", nrow(mdata), "subjects!", sep = " "))
  }

  ##random effect covariates
  if (model == "interslope") {
    if (prod(RE %in% ynames) == 0) {
      Fakename <- which(RE %in% ynames == FALSE)
      stop(paste0("The variable ", RE[Fakename], " not found"))
    } else {
      p1a <- 1 + length(RE)
      random.xnam <- paste(RE[1:length(RE)], sep = "")
      fmla.random <- paste("(", paste(random.xnam, collapse= "+"), "|", ID, ")", sep = "")
      ydata_RE <- data.frame(1, ydata[, RE])
      colnames(ydata_RE)[2:length(ydata_RE)] <- paste0(RE, "_RE")
      colnames(ydata_RE)[1] <- "intercept_RE"
    }
  } else if (model == "intercept") {
    if (!is.null(RE)) {
      stop("You are fitting a mixed effects model with random intercept only
           but random effects covariates are specified at the same time. Please respecify your model!")
    }
    fmla.random <- paste("(", 1, "|", ID, ")", sep = "")
    p1a <- 1
    ydata_RE <- as.data.frame(rep(1, n1))
    colnames(ydata_RE) <- "intercept_RE"
  } else {
    stop("model should be one of the following options: interslope or intercept.")
  }

  ##fixed effects covariates
  ydata_FE <- data.frame(1, ydata[, long[2:length(long)]])
  p1 <- ncol(ydata_FE)
  colnames(ydata_FE)[1] <- "intercept"
  colnames(ydata_FE)[2:ncol(ydata_FE)] <- long[2:length(long)]
  xnam <- paste(long[2:length(long)], sep = "")
  fmla.fixed <- paste(long[1], " ~ ", paste(xnam, collapse= "+"))
  fmla <- as.formula(paste(fmla.fixed, fmla.random, sep = "+"))

  writeLines("Model fit is starting!")
  lmem <- lme4::lmer(formula = fmla, data = ydata, REML = REML)
  yfile <- cbind(ydata[, long[1]], ydata_RE, ydata_FE)
  colnames(yfile)[1] <- long[1]

  cfile <- cdata[, survival]

  if (sum(complete.cases(yfile)) != n1) {
    stop("Missing values detected in ydata! Please make sure your longitudinal data is complete!")
  }

  if (sum(complete.cases(cfile)) != cdim[1]) {
    stop("Missing values detected in cdata! Please make sure your survival data is complete!")
  }

  if (p1a > 3) {
    stop("Maximum of 3 random effects are allowed. Please reconsider the random effect covariates you need!")
  }

  if((p1<1)|(p1a<1)){
    stop("Possibe wrong dimension of fixed effects in Y!")
  }

  # number of subjects in study is equals to the #of rows in C matrix
  k=cdim[1]
  p2=length(survival)-2

  yfilenew=tempfile(pattern = "", fileext = ".txt")
  writenh(yfile,yfilenew)
  mfilenew=tempfile(pattern = "", fileext = ".txt")
  writenh(mdata,mfilenew)
  cfilenew=tempfile(pattern = "", fileext = ".txt")
  writenh(cfile,cfilenew)

  D <- lme4::VarCorr(lmem)
  name <- names(D)
  D <- as.data.frame(D[name])
  D <- as.matrix(D)
  Sigcovfile = tempfile(pattern = "", fileext = ".txt")
  writenh(D, Sigcovfile)
  ##Residuals
  sigma <- sigma(lmem)^2
  ##marginal covariance matrix V
  beta <- lme4::fixef(lmem)
  beta <- c(beta, sigma)
  #beta <- c(91.314884, 0.350534, 0.706091, 4.157956, -2.360598, 2.630575, -2.584545, -1.674665, 168.939393)
  beta <- matrix(beta, ncol = 1, nrow = length(beta))
  Betasigmafile = tempfile(pattern = "", fileext = ".txt")
  writenh(beta, Betasigmafile)

  if (prod(c(1, 2) %in% unlist(unique(cfile[, 2]))) == 1) {
    if (survinital) {
      ## fit a Cox model
      surv_xnam <- paste(survival[3:length(survival)], sep = "")
      survfmla.fixed <- paste(surv_xnam, collapse= "+")
      survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==1)")
      survfmla <- as.formula(paste(survfmla.out1, survfmla.fixed, sep = "~"))
      fitSURV1 <- survival::coxph(formula = survfmla, data = cfile, x = TRUE)
      gamma1 <- fitSURV1$coefficients

      survfmla.out2 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==2)")
      survfmla <- as.formula(paste(survfmla.out2, survfmla.fixed, sep = "~"))
      fitSURV2 <- survival::coxph(formula = survfmla, data = cfile, x = TRUE)
      gamma2 <- fitSURV2$coefficients
    } else {
      gamma1 <- rep(0, p2)
      gamma2 <- rep(0, p2)
    }

    myresult=jmcs_main(tol, k, n1, p1, p2, p1a, maxiter, point, xs, ws, yfilenew, cfilenew,
                       mfilenew, Betasigmafile, Sigcovfile, gamma1, gamma2, trace)
    myresult$type="jmcs";
  } else if (prod(c(1) %in% unique(cfile[, 2])) == 1) {

    if (survinital) {
      ## fit a Cox model
      surv_xnam <- paste(survival[3:length(survival)], sep = "")
      survfmla.fixed <- paste(surv_xnam, collapse= "+")
      survfmla.out <- paste0("survival::Surv(", survival[1], ", ", survival[2], ")")
      survfmla <- as.formula(paste(survfmla.out, survfmla.fixed, sep = "~"))
      fitSURV <- survival::coxph(formula = survfmla, data = cfile, x = TRUE)
      gamma <- fitSURV$coefficients
    } else {
      gamma <- rep(0, p2)
    }

    myresult=jmcsf_main(tol, k, n1, p1, p2, p1a, maxiter, point, xs, ws, yfilenew, cfilenew,
                        mfilenew, Betasigmafile, Sigcovfile, gamma, trace)
    myresult$type="jmcsf";
  } else {
    stop(paste0(unique(cfile[, 2]), " is not an appropriate code of single / competing risks failure type.
                  Please correctly specify the event variable."))
  }
  PropComp <- round(table(cfile[, 2])/k * 100, 2)
  ynames=colnames(yfile)
  ydim=dim(yfile)
  cnames=colnames(cfile)
  #names
  if (!is.null(myresult$betas)) {
    names(myresult$betas)=ynames[(1+p1a+1):ydim[2]]
  }
  if (!is.null(myresult$gamma_matrix)) {
    colnames(myresult$gamma_matrix)=cnames[3:(p2+3-1)]
  }

  myresult$k=k

  LongOut <- ynames[1]
  LongX <- paste0(ynames[(3+p1a):length(ynames)], collapse = "+")
  FunCall_long <- as.formula(paste(LongOut, LongX, sep = "~"))

  ##create survival submodel formula
  #cnames <- colnames(cdata)
  SurvOut <- paste0("Surv(", cnames[1], ",", cnames[2], ")")
  SurvX <- paste0(cnames[-(1:2)], collapse = "+")
  FunCall_survival <- as.formula(paste(SurvOut, SurvX, sep = "~"))


  DataPath <- NULL
  SummaryInfo <- list(k, n1, PropComp, FunCall_long, FunCall_survival, DataPath)
  names(SummaryInfo) <- c("NumSub", "Numobs", "PropComp",
                          "LongitudinalSubmodel", "SurvivalSubmodel", "DataPath")

  myresult$SummaryInfo <- SummaryInfo
  myresult$point <- point
  myresult$REML <- REML
  myresult$ydata <- ydata
  myresult$cdata <- cdata
  myresult$mdata <- mdata

  class(myresult) <- "FastJM"

  return (myresult)

}
