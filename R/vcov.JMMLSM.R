##' @title Variance-covariance matrix of the estimated parameters for joint models
##' @name vcov
##' @aliases vcov.JMMLSM
##' @description Extract variance-covariance matrix for joint models.
##' @param object an object inheriting from class \code{JMMLSM}.
##' @param ... further arguments passed to or from other methods.
##' @return a matrix of variance covariance of all parameter estimates.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{JMMLSM}}
##' @export
##' 
vcov.JMMLSM <- function(object, ...) {
  if (!inherits(object, "JMMLSM"))
    stop("Use only with 'JMMLSM' objects.\n")
  
  # ----- helpers -------------------------------------------------------------
  make_sig_names_ordered <- function(p, prefix = "Sig") {
    # Diagonals first: Sig11, Sig22, ..., Sigpp
    out <- sprintf("%s%d%d", prefix, seq_len(p), seq_len(p))
    # Then off-diagonals by increasing distance d = j - i
    for (d in 1:(p - 1)) {
      for (i in 1:(p - d)) {
        j <- i + d
        out <- c(out, sprintf("%s%d%d", prefix, i, j))
      }
    }
    out
  }
  
  make_alpha_labels <- function(p1a, random, suffix) {
    base <- c("(Intercept)", random)
    if (length(base) < (p1a - 1)) {
      base <- c(base, paste0("RE", seq_len((p1a - 1) - length(base))))
    }
    base <- base[seq_len(p1a - 1)]
    paste0("T.asso:", base, suffix)
  }
  # --------------------------------------------------------------------------
  
  variance.formula <- as.formula(paste("", object$LongitudinalSubmodelvariance[3], sep = "~"))
  
  getdum <- getdummy.JMH(
    long.formula = object$LongitudinalSubmodelmean,
    surv.formula = object$SurvivalSubmodel,
    variance.formula = variance.formula,
    random = object$random, ydata = object$ydata, cdata = object$cdata
  )
  
  random <- all.vars(object$random)
  vcov   <- as.data.frame(object$vcov)
  
  long <- paste0("Ymean.", names(object$beta))
  tau  <- paste0("Yvariance.", names(object$tau))
  
  survival1 <- paste0("T.", names(object$gamma1))
  
  # p1a = number of association terms + intercept
  p1a <- length(object$alpha1) + 1L
  
  # Association labels (event 1)
  alpha1 <- make_alpha_labels(p1a, random, "_1")
  nu1    <- "T.asso.Var_(Intercept)_1:"
  
  # Sig names in your requested order
  sig <- make_sig_names_ordered(p1a, prefix = "Sig")
  
  if (isTRUE(object$CompetingRisk)) {
    survival2 <- paste0("T.", names(object$gamma2))
    alpha2    <- make_alpha_labels(p1a, random, "_2")
    nu2       <- "T.asso.Var_(Intercept)_2:"
    
    colnames(vcov) <- c(long, tau, survival1, survival2, alpha1, alpha2, nu1, nu2, sig)
    rownames(vcov) <- c(long, tau, survival1, survival2, alpha1, alpha2, nu1, nu2, sig)
  } else {
    colnames(vcov) <- c(long, tau, survival1, alpha1, nu1, sig)
    rownames(vcov) <- c(long, tau, survival1, alpha1, nu1, sig)
  }
  
  vcov
}