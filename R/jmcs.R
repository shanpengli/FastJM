##' Joint modeling of longitudinal continuous data and competing risks
##' @title Joint Modelling for Continuous outcomes
##' @param ydata a longitudinal data frame in long format.
##' @param cdata a survival data frame with competing risks or single failure.
##' Each subject has one data entry.
##' @param long.formula a formula object with the response variable and fixed effects covariates
##' to be included in the longitudinal sub-model.
##' @param surv.formula a formula object with the survival time, event indicator, and the covariates
##' to be included in the survival sub-model.
##' @param ID an element with the name identification of subject to
##' identify subject identification in the dataframes.
##' @param RE a vector of the names of random effect covariates to be included
##' in the longitudinal sub-model.
##' @param model an element that specifies the model type of random effects.
##' Deafult is "interslope", i.e., both random intercept and slope are included.
##' @param point the number of pseudo-adaptive Gauss-Hermite quadrature points
##' to be chosen for numerical integration. Default is 6 which produces stable estimates in most dataframes.
##' @param maxiter the maximum number of iterations of the EM algorithm that the function will perform. Default is 10000.
##' @param do.trace Print detailed information of each iteration. Default is FALSE, i.e., not to print the iteration details.
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

jmcs <- function(ydata, cdata, long.formula, surv.formula, ID, RE, model = "interslope",
                 point = 6, maxiter = 10000, do.trace = FALSE)
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
  if (is.null(model)) {
    stop("model type must be specified. It should be one of the following options: intercept or interslope.")
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
  xnam <- paste(long[2:length(long)], sep = "")
  fmla.fixed <- paste(long[1], " ~ ", paste(xnam, collapse= "+"))
  fmla <- as.formula(paste(fmla.fixed, fmla.random, sep = "+"))

  writeLines("Model fit is starting!")
  lmem <- lme4::lmer(formula = fmla, data = ydata, REML = FALSE)
  yfile <- cbind(ydata[, long[1]], ydata_RE, ydata_FE)
  colnames(yfile)[1] <- long[1]

  cfile <- cdata[, survival]

  if (sum(complete.cases(yfile)) != n1) {
    stop("Missing values detected! Please make sure your longitudinal data is complete!")
  }

  if (sum(complete.cases(cfile)) != cdim[1]) {
    stop("Missing values detected! Please make sure your survival data is complete!")
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
  beta <- matrix(beta, ncol = 1, nrow = length(beta))
  Betasigmafile = tempfile(pattern = "", fileext = ".txt")
  writenh(beta, Betasigmafile)
  if (prod(c(0, 1, 2) %in% unique(cfile[, 2])) == 1) {
    myresult=jmcs_main(k, n1, p1, p2, p1a, maxiter, point, xs, ws, yfilenew, cfilenew,
                       mfilenew, Betasigmafile, Sigcovfile, trace)
    myresult$type="jmcs";
  } else if (prod(c(0, 1) %in% unique(cfile[, 2])) == 1) {
    myresult=jmcsf_main(k, n1, p1, p2, p1a, maxiter, point, xs, ws, yfilenew, cfilenew,
                        mfilenew, Betasigmafile, Sigcovfile, trace)
    myresult$type="jmcsf";
  } else {
    stop(paste0(unique(cfile[, 2]), " is not an appropriate code of single / competing risks failure type.
                  Please correctly specify the event variable."))
  }
  ynames=colnames(yfile)
  ydim=dim(yfile)
  cnames=colnames(cfile)
  #names
  names(myresult$betas)=ynames[(1+p1a+1):ydim[2]]

  colnames(myresult$gamma_matrix)=cnames[3:(p2+3-1)]

  myresult$k=k
  class(myresult) <- "FastJM"

  return (myresult)

}
