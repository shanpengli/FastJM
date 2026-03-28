##' @title Plot conditional probabilities for new subjects
##' @name plot.survfitJMMLSM
##' @aliases plot.survfitJMMLSM
##' @description Plot conditional probabilities for new subjects. 
##' If \code{CompetingRisk = FALSE}, print the survival probabilities. 
##' Otherwise, print the cumulative incidence probabilities for each failure type.
##' @param x x of class \code{survfitJMMLSM}.
##' @param include.y include longitudinal responses of this subject versus time. Default is FALSE.
##' @param xlab X axis label.
##' @param ylab Y axis label.
##' @param xlim X axis support.
##' @param ylim.long Y axis support for the longitudinal outcome.
##' @param ylim.surv Y axis support for the event / survival probability.
##' @param ... further arguments passed to or from other methods.
##' @return plots of conditional probabilities over different pre-specified time points for subjects. 
##' If single failure type, then survival probabilities will be returned. 
##' Otherwise, cumulative incidence probabilities for each failure type will be returned.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{survfitJMMLSM}}
##' @export
##' 
plot.survfitJMMLSM <- function (x, include.y = FALSE, xlab = NULL, ylab = NULL, 
                              xlim = NULL, ylim.long = NULL, ylim.surv = NULL, ...) {
  
  if (!inherits(x, "survfitJMMLSM"))
    stop("Use only with 'survfitJMMLSM' xs.\n")
  
  if (is.null(xlab)) xlab = "Time"
  if (is.null(ylim.surv)) ylim.surv = c(0, 1)
  
  if (!x$CompetingRisk) {
    which = 1:nrow(x$Last.time)
    ask = (prod(par("mfcol")) < length(which))
    show <- rep(TRUE, length(which))
    return = FALSE
    
    if (!return) {
      one.fig <- prod(par("mfcol")) == 1 
    }
    if (ask && !return) {
      op <- par(ask = TRUE)
      on.exit(par(op))
    }
    
    for (i in 1:nrow(x$Last.time)) {
      if (show[i] && !return) {
        times <- c(0, as.numeric(x$Last.time[i, 2]), x$Pred[[i]][, 1])
        probmean <- c(1, 1, x$Pred[[i]][, 2])
        if (is.null(xlim)) xlim <- c(0, max(x$Pred[[i]][, 1]))
        
        if (!include.y) {
          if (is.null(ylab)) ylab <- expression(paste("Pr(", T[i] >= u, " | ", T[i] > s, 
                                                      ", ", y[i]^(s), ", ",  Psi,")", sep = " "))
          plot(times, probmean, xlab = xlab, ylab = ylab, xlim = xlim,
               main = paste("Subject", x$Last.time[i, 1], sep = " "), col = "red", type = "l", ylim = ylim.surv)
          segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                   y1 = 1,
                   lwd = 1)
        } else {
          if (is.null(ylab)) ylab = "Longitudinal outcome"
          plot(x$y.obs[[i]][, 1], x$y.obs[[i]][, 2], xlim = xlim, axes = TRUE, xlab = xlab, 
               ylab = "", type = "p", pch = 8, ylim = ylim.long)
          title(ylab = ylab, line=2.5)
          par(new = TRUE)    
          plot(times, probmean, xlab = "", ylab = "", 
               main = paste("Subject", x$Last.time[i, 1], sep = " "), xlim = xlim, col = "red", type = "l", 
               ylim = ylim.surv, axes = FALSE)
          axis(side = 4, at = pretty(range(ylim.surv)), line = 0) 
          mtext(expression(paste("Pr(", T[i] >= u, " | ", T[i] > s, 
                                 ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), side = 4, line = 2.5)
          segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                   y1 = 1,
                   lwd = 1)
        }
      
      }
    }
    invisible()
    
  } else {
    
    which = 1:nrow(x$Last.time)
    ylim <- c(0, 1)
    
    ask = (prod(par("mfcol")) < length(which))
    show <- rep(TRUE, 2*length(which))
    return = FALSE
    
    if (!return) {
      one.fig <- prod(par("mfcol")) == 1 
    }
    if (ask && !return) {
      op <- par(ask = TRUE)
      on.exit(par(op))
    }
    
    for (i in 1:nrow(x$Last.time)) {
      for (j in 1:2) {
        if (is.null(xlim)) xlim <- c(0, max(x$Pred[[i]][, 1]))
        if (show[(i-1)*2+j] && !return) {
          times <- c(0, as.numeric(x$Last.time[i, 2]), x$Pred[[i]][, 1])
          probmean <- c(0, 0, x$Pred[[i]][, j+1])
          
          if (!include.y) {
            plot(times, probmean, xlab = xlab, ylab = expression(paste("Pr(", T[i] <= u, ",", D[i] == k, " | ", T[i] > s, 
                                                                       ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), 
                 main = paste("Subject", x$Last.time[i, 1], "k =", j, sep = " "), 
                 col = "red", type = "l", ylim = ylim.surv, xlim = xlim)
            segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                     y1 = 1,
                     lwd = 1)
          } else {
            if (is.null(ylab)) ylab = "Longitudinal outcome"
            plot(x$y.obs[[i]][, 1], x$y.obs[[i]][, 2], xlim = xlim, axes = TRUE, 
                 xlab = xlab, ylab = "", type = "p", pch = 8, ylim = ylim.long)
            title(ylab = ylab, line=2.5)
            par(new = TRUE)    
            plot(times, probmean, xlab = "", ylab = "", 
                 main = paste("Subject", x$Last.time[i, 1], "k =", j, sep = " "), 
                 col = "red", type = "l", ylim = ylim.surv, axes = FALSE, xlim = xlim)
            axis(side = 4, at = pretty(range(ylim.surv)), line = 0) 
            mtext(expression(paste("Pr(", T[i] <= u, ",", D[i] == k, " | ", 
                                   T[i] > s, ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), 
                  side = 4, line = 2.5)
            segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                     y1 = 1,
                     lwd = 1)
          }
        
        }
      }
    }
    invisible()
  }
  
}