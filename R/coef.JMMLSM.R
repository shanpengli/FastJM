##' @rdname coef
##' @export

coef.JMMLSM <- function(object, ...) {
  if (!inherits(object, "JMMLSM")) stop("Use only with 'JMMLSM' objects.")
  
  named_block <- function(x, nm) {
    if (length(x) != length(nm) || anyNA(nm) || any(!nzchar(nm)))
      stop("Missing names or inconsistent parameter dimensions.")
    stats::setNames(as.numeric(x), nm)
  }
  coefficient <- function(x, prefix) {
    if (is.null(names(x))) stop("Regression coefficient names are missing.")
    named_block(x, paste0(prefix, names(x)))
  }
  if (!is.logical(object$CompetingRisk) || length(object$CompetingRisk) != 1L ||
      is.na(object$CompetingRisk))
    stop("Invalid 'CompetingRisk' flag.")
  cr <- object$CompetingRisk
  
  required <- c("beta", "tau", "gamma1", "alpha1", "vee1", "Sig")
  if (cr) required <- c(required, "gamma2", "alpha2", "vee2")
  
  valid <- vapply(required, function(nm) {
    x <- object[[nm]]
    is.numeric(x) && length(x) > 0L && all(is.finite(x))
  }, logical(1))
  if (!all(valid)) stop("Missing or invalid parameter estimates: ",
                        paste(required[!valid], collapse = ", "), ".")
  S <- as.matrix(object$Sig)
  q <- nrow(S)
  if (q < 1L || ncol(S) != q) stop("'Sig' must be a nonempty square matrix.")
  
  n_mean <- length(object$alpha1)
  if (q != n_mean + 1L || length(object$vee1) != 1L ||
      (cr && (length(object$alpha2) != n_mean || length(object$vee2) != 1L)))
    stop("Inconsistent association or random-effects dimensions.")
  re <- c("(Intercept)", all.vars(object$random))[seq_len(n_mean)]
  if (anyNA(re)) stop("Cannot identify random-effect names.")
  ans <- c(coefficient(object$beta, "Ymean."),
           coefficient(object$tau, "Yvariance."),
           coefficient(object$gamma1, "T."))
  if (cr) ans <- c(ans, coefficient(object$gamma2, "T."))
  ans <- c(ans, named_block(object$alpha1, paste0("T.asso:", re, "_1")))
  if (cr) ans <- c(ans, named_block(object$alpha2, paste0("T.asso:", re, "_2")))
  ans <- c(ans, named_block(object$vee1, "T.asso.Var_(Intercept)_1:"))
  if (cr) ans <- c(ans, named_block(object$vee2, "T.asso.Var_(Intercept)_2:"))
  cov_names <- function(i, j) paste0("Sig", i, j)
  
  covariance <- named_block(diag(S), cov_names(seq_len(q), seq_len(q)))
  for (lag in seq_len(q - 1L)) {
    i <- seq_len(q - lag)
    j <- i + lag
    covariance <- c(covariance, named_block(S[cbind(i, j)], cov_names(i, j)))
  }
  c(ans, covariance)
}
