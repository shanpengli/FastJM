##' @rdname coef
##' @export

coef.mvjmcs <- function(object, ...) {
  if (!inherits(object, "mvjmcs")) stop("Use only with 'mvjmcs' objects.")
  
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
  
  required <- c("beta", "sigma", "gamma1", "alpha1", "Sig")
  if (cr) required <- c(required, "gamma2", "alpha2")
  
  valid <- vapply(required, function(nm) {
    x <- object[[nm]]
    is.numeric(x) && length(x) > 0L && all(is.finite(x))
  }, logical(1))
  if (!all(valid)) stop("Missing or invalid parameter estimates: ",
                        paste(required[!valid], collapse = ", "), ".")
  S <- as.matrix(object$Sig)
  q <- nrow(S)
  if (q < 1L || ncol(S) != q) stop("'Sig' must be a nonempty square matrix.")
  
  n_bio <- length(object$sigma)
  if (!is.list(object$random) || length(object$random) != n_bio)
    stop("Expected one random-effects formula per biomarker.")
  re <- cov_re <- character(0)
  for (g in seq_len(n_bio)) {
    vars <- all.vars(object$random[[g]])
    if (!length(vars)) stop("Invalid random-effects formula.")
    vars <- vars[-length(vars)]  # Final variable is the grouping variable.
    re <- c(re, paste0(c("(Intercept)", vars), "_bio", g))
    cov_re <- c(cov_re, paste0(c("Intercept", vars), g))
  }
  if (length(cov_re) != q)
    stop("Random-effects formulas do not match 'Sig'.")
  association <- object$latAsso
  if (is.null(association)) association <- "sre"
  if (length(association) != 1L || is.na(association) ||
      !association %in% c("sre", "present", "presentlp"))
    stop("Unknown latent association structure.")
  assoc_names <- if (association == "sre") re else paste0(association, "_bio", seq_len(n_bio))
  ans <- c(coefficient(object$beta, "Y."),
           named_block(object$sigma, paste0("Y.sigma^2_bio", seq_len(n_bio))),
           coefficient(object$gamma1, "T."))
  if (cr) ans <- c(ans, coefficient(object$gamma2, "T."))
  ans <- c(ans, named_block(object$alpha1, paste0("T.asso:", assoc_names, "_1")))
  if (cr) ans <- c(ans, named_block(object$alpha2, paste0("T.asso:", assoc_names, "_2")))
  cov_names <- function(i, j) paste0(cov_re[i], ":", cov_re[j])
  
  covariance <- named_block(diag(S), cov_names(seq_len(q), seq_len(q)))
  for (lag in seq_len(q - 1L)) {
    i <- seq_len(q - lag)
    j <- i + lag
    covariance <- c(covariance, named_block(S[cbind(i, j)], cov_names(i, j)))
  }
  c(ans, covariance)
}
