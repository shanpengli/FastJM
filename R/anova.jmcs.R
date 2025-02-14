##' @title Anova Method for Fitted Joint Models
##' @name anova
##' @aliases anova.jmcs
##' @description Performs a likelihood ratio test between two nested joint models.
##' @param object an object inheriting from class \code{jmcs}, nested in \code{object2}.
##' @param object2 an object inheriting from class \code{jmcs}.
##' @param digits the number of significant digits to use when printing. Default is 4.
##' @param ... further arguments passed to or from other methods.
##' @return A table to summarize the likelihood ratio test.
##' @seealso \code{\link{jmcs}}
##' @examples 
##' \donttest{
##' # Fit a joint model
##' fit <- jmcs(ydata = ydata, cdata = cdata, 
##'             long.formula = response ~ time + x1, 
##'             surv.formula = Surv(surv, failure_type) ~ x1 + x2, 
##'             random =  ~ time| ID)
##'
##'fit2 <- jmcs(ydata = ydata, cdata = cdata, 
##'             long.formula = response ~ time + gender + x1 + race, 
##'             surv.formula = Surv(surv, failure_type) ~ x1 + gender + x2 + race, 
##'             random =  ~ time| ID)      
##'
##'anova(fit, fit2)   
##'}
##'              
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @export
##' 

anova.jmcs <- function(object, object2, digits = 4, ...) {
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  
  if (!inherits(object2, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  
  if (object$CompetingRisk != object2$CompetingRisk)
    stop("The survival submodels of these two objects are differently structured and thus not available for anova.")
  
  ### data log likelihood
  L0 <- object$loglike
  L1 <- object2$loglike
  
  if (object$CompetingRisk) {
    df0 <- length(object$beta) + 2*(length(object$gamma1) + length(object$nu1)) + 
      ncol(object$Sig)*(ncol(object$Sig) + 1)/2 + 1
    df1 <- length(object2$beta) + 2*(length(object2$gamma1) + length(object2$nu1)) + 
      ncol(object2$Sig)*(ncol(object2$Sig) + 1)/2 + 1
  } else {
    df0 <- length(object$beta) + length(object$gamma1) + length(object$nu1) + 
      ncol(object$Sig)*(ncol(object$Sig) + 1)/2 + 1
    df1 <- length(object2$beta) + length(object2$gamma1) + length(object2$nu1) + 
      ncol(object2$Sig)*(ncol(object2$Sig) + 1)/2 + 1
  }
  
  AIC0 <- -2*L0 + 2*df0
  AIC1 <- -2*L1 + 2*df1
  BIC0 <- -2*L0 + log(nrow(object$cdata))*df0
  BIC1 <- -2*L1 + log(nrow(object$cdata))*df1
  
  dff <- df1 - df0
  
  if (dff < 0) {
    stop("'object' should be nested in 'object2'.")
  }
  
  LRT <- -2*(L0 - L1)
  
  if (LRT < 0) {
    stop("'object' is either not nested in 'object2' or 'object2' has numerical issue.")
  }
  
  pval <- pchisq(LRT, dff, lower.tail = FALSE)
  
  table <- matrix(" ", nrow = 2, ncol = 6)
  table <- as.data.frame(table)
  rownames(table)[1] <- deparse(substitute(object))
  rownames(table)[2] <- deparse(substitute(object2))
  colnames(table) <- c("AIC", "BIC", "loglike", "LRT", "df", "p-value")
  table[1, 1:3] <- c(AIC0, BIC0, L0)
  table[1, 1:3] <- sprintf(paste("%.", digits, "f", sep = ""), table[1, 1:3])
  table[2, 1:6] <- c(AIC1, BIC1, L1, LRT, dff, pval)
  table[2, c(1:4, 6)] <- sprintf(paste("%.", digits, "f", sep = ""), table[2, c(1:4, 6)])
  table
}