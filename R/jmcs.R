##' Joint modeling of longitudinal continuous data and competing risks
##' @title Joint Modelling for Continuous outcomes
##' @param p1  The dimension of fixed effects (including intercept) in yfile.
##' @param yfile Y matrix for longitudinal measurements in long format. For example, for a subject with n measurements, there should be n rows for this subject. The # of rows in y matrix is the total number of measurements for all subjects in the study. The columns in Y should start with the longitudinal outcome (column 1), the covariates for the random effects, and then the covariates for the fixed effects.
##' @param cfile C matrix for competing risks failure time data. Each subject has one data entry, so the number of rows equals to the number of subjects. The survival / censoring time is included in the first column, and the failure type coded as 0 (censored events), 1 (risk 1), or 2 (risk 2) is given in the second column. Two competing risks are assumed. The covariates are included in the third column and on.
##' @param mfile M vector to indicate the number of longitudinal measurements per subject. The number of rows equals to the number of subjects.
##' @param point Quadrature points used in the EM procedure.Default is 20.
##' @param maxiter Maximum values of iterations. Default is 100000.
##' @param do.trace Print detailed information of each iteration. Default is false, i.e., not to print the iteration details.
##' @param type_file Types of inputs. Default is true, i.e.  data files with headers. If set to "F", inputs are changed to data matrixes or data.frames (with headers)
##' @param ... further arguments passed to or from other methods.
##' @return Object of class \code{FastJM} with elements
##'   \tabular{ll}{
##'       \code{vcmatrix}    \tab  The variance-covariance matrix for all the parameters. The parameters are in the order: \eqn{\beta}, \eqn{\sigma^2}, \eqn{\gamma}, \eqn{\nu}, and \eqn{\Sigma}. The elements in \eqn{\Sigma} are output in the order along the main diagonal line, then the second main diagonal line, and so on. \cr
##'       \code{betas} \tab The point  estimates of \eqn{\beta}. \cr
##'       \code{se_betas} \tab The standard error estimate of \eqn{\beta}. \cr
##'       \code{gamma_matrix} \tab  The point  estimate of \eqn{\gamma}. \cr
##'       \code{se_gamma_matrix}   \tab  The standard error estimate of \eqn{\gamma}. \cr
##'       \code{v_estimate} \tab The point  estimate of \eqn{\nu}. \cr
##'       \code{se_v_estimate}    \tab The standard error estimate of \eqn{\nu}. \cr
##'       \code{sigma2_val}     \tab  The point estimate of \eqn{\sigma^2}.\cr
##'       \code{se_sigma2_val}     \tab  The standard error estimate of \eqn{\sigma^2}.\cr
##'       \code{sigma_matrix}     \tab The point estimate of \eqn{\Sigma} (only the upper triangle portion of the matrix is output).\cr
##'       \code{se_sigma}     \tab The standard error estimate of \eqn{\Sigma}.The standard errors are given in this order: main diagonal, the second main diagonal, and so on. \cr
##'       \code{loglike}     \tab Log Likelihood.\cr
##'   }
##'
##' @examples
##' # A toy example on simulated data
##' require(FastJM)
##' set.seed(123)
##' yfile=system.file("extdata", "simy0.txt", package = "FastJM")
##' cfile=system.file("extdata", "simc0.txt", package = "FastJM")
##' mfile=system.file("extdata", "simm0.txt", package = "FastJM")
##' res2=jmcs(p1=3,yfile,cfile,mfile,point=6,type_file=TRUE)

##' @references
##' \itemize{
##' \item Elashoff, Robert M., Gang Li, and Ning Li. "A joint model for longitudinal measurements and survival data in the presence of multiple failure types." Biometrics 64.3 (2008): 762-771.
##' }
##' @export

jmcs<- function (p1,yfile,cfile,mfile,point=6,maxiter=10000,do.trace=FALSE,type_file=TRUE)
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

  if (type_file){
    # store header names for future useage

    ydata=read.table(yfile)
    ynames=colnames(ydata)

    cdata=read.table(cfile)
    cnames=colnames(cdata)
    cfile=tempfile(pattern = "", fileext = ".txt")
    writenh(cdata,cfile)
    mdata=read.table(mfile)
    mfilenew=tempfile(pattern = "", fileext = ".txt")
    writenh(mdata,mfilenew)

    ydim=dim(ydata)
    # number of observations in study is equals to the #of rows in Y matrix
    n1=ydim[1]
    # number of random effects
    p1a=ydim[2]-1-p1-1

    if((p1<1)|(p1a<1)){
      stop("Possibe wrong dimension of fixed effects in Y!")
    }

    cdim=dim(cdata);

    # number of subjects in study is equals to the #of rows in C matrix
    k=cdim[1]
    p2=cdim[2]-2
    colnames(ydata)[1] <- "ID"

    ## fit a linear mixed effect model
    name <- colnames(ydata)[(4+p1a):ncol(ydata)]
    xnam <- paste(name[1:length(name)], sep = "")
    random.name <- colnames(ydata)[4:(2+p1a)]
    random.xnam <- paste(random.name[1:length(random.name)], sep = "")

    fmla.fixed <- paste(colnames(ydata)[2], " ~ ", paste(xnam, collapse= "+"))
    fmla.random <- paste("(", paste(random.xnam, collapse= "+"), "|", colnames(ydata)[1], ")", sep = "")
    fmla <- as.formula(paste(fmla.fixed, fmla.random, sep = "+"))
    writeLines("Model fit is starting!")
    lmem <- lme4::lmer(formula = fmla, data = ydata, REML = FALSE)

    ydata <- ydata[, -1]
    yfile=tempfile(pattern = "", fileext = ".txt")
    writenh(ydata,yfile)

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

    if (prod(c(0, 1, 2) %in% unique(cdata[, 2])) == 1) {
      myresult=jmcs_main(k, n1, p1, p2, p1a, maxiter, point, xs, ws, yfile, cfile,
                         mfilenew, Betasigmafile, Sigcovfile, trace)
      myresult$type="jmcs";
    } else if (prod(c(0, 1) %in% unique(cdata[, 2])) == 1) {
      myresult=jmcsf_main(k, n1, p1, p2, p1a, maxiter, point, xs, ws, yfile, cfile,
                         mfilenew, Betasigmafile, Sigcovfile, trace)
      myresult$type="jmcsf";
    } else {
      stop(paste0(unique(cdata[, 2]), " is not an appropriate code of single / competing risks failure type.
                  Please correctly specify the event variable."))
    }


  }else{

    cnames=colnames(cfile)
    cfilenew=tempfile(pattern = "", fileext = ".txt")
    writenh(cfile,cfilenew)

    mfilenew=tempfile(pattern = "", fileext = ".txt")
    writenh(mfile,mfilenew)

    ydim=dim(yfile)
    # number of observations in study is equals to the #of rows in Y matrix
    n1=ydim[1]
    # number of random effects
    p1a=ydim[2]-1-p1-1
    if((p1<1)|(p1a<1)){
      stop("Possibe wrong dimension of fixed effects in Y!")
    }

    cdim=dim(cfile)

    # number of subjects in study is equals to the #of rows in C matrix
    k=cdim[1]
    p2=cdim[2]-2
    colnames(yfile)[1] <- "ID"

    ## fit a linear mixed effect model
    name <- colnames(yfile)[(4+p1a):ncol(yfile)]
    xnam <- paste(name[1:length(name)], sep = "")
    random.name <- colnames(yfile)[4:(2+p1a)]
    random.xnam <- paste(random.name[1:length(random.name)], sep = "")

    fmla.fixed <- paste(colnames(yfile)[2], " ~ ", paste(xnam, collapse= "+"))
    fmla.random <- paste("(", paste(random.xnam, collapse= "+"), "|", colnames(yfile)[1], ")", sep = "")
    fmla <- as.formula(paste(fmla.fixed, fmla.random, sep = "+"))
    writeLines("Model fit is starting!")
    lmem <- lmer(formula = fmla, data = yfile, REML = FALSE)

    yfile <- yfile[, -1]
    yfilenew=tempfile(pattern = "", fileext = ".txt")
    writenh(yfile,yfilenew)
    D <- VarCorr(lmem)
    name <- names(D)
    D <- as.data.frame(D[name])
    D <- as.matrix(D)
    Sigcovfile = tempfile(pattern = "", fileext = ".txt")
    writenh(D, Sigcovfile)
    ##Residuals
    sigma <- sigma(lmem)^2
    ##marginal covariance matrix V
    beta <- fixef(lmem)
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

  }


  #names
  ynames=colnames(yfile)
  ydim=dim(yfile)
  names(myresult$betas)=ynames[(1+p1a+1):ydim[2]]

  colnames(myresult$gamma_matrix)=cnames[3:(p2+3-1)]

  myresult$k=k
  class(myresult) <- "FastJM"

  return (myresult)
}
