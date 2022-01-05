##' @title Plot survival probabilities for new subjects
##' @name plot.survfitjmcs
##' @aliases plot.survfitjmcs
##' @description Plot survival probabilities for new subjects.
##' @param x x of class \code{survfitjmcs}.
##' @param estimator character string specifying, whether to include in the plot 
##' the mean of the conditional probabilities of survival, the median or both. 
##' The mean and median are taken as estimates of these conditional probabilities 
##' over the M replications of the Monte Carlo scheme described in \code{survfitjmcs}.
##' @param conf.int logical; if TRUE, then a pointwise confidence interval is included in the plot. Default is TRUE.
##' @param include.y include longitudinal responses of this subject versus time. Default is FALSE.
##' @param ... further arguments passed to or from other methods.
##' @return plots of conditional probabilities over different pre-specified time points for subjects. 
##' If single failure type, then survival probabilities will be returned. 
##' Otherwise, cumulative incidence probabilities for each failure type will be returned.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{survfitjmcs}}
##' @examples 
##' \donttest{
##' # a joint model fit
##' fit <- jmcs(ydata = ydata, cdata = cdata, 
##' long.formula = response ~ time + x1, 
##' surv.formula = Surv(surv, failure_type) ~ x1 + x2, 
##' random =  ~ time| ID)
##' 
##' # Prediction of cumulative incidence for competing risks data
##' ND <- ydata[ydata$ID %in% c(419, 218), ]
##' ID <- unique(ND$ID)
##' NDc <- cdata[cdata$ID  %in% ID, ]
##' survfit <- survfitjmcs(fit, 
##'                        ynewdata = ND, 
##'                        cnewdata = NDc, 
##'                        u = seq(3, 4.8, by = 0.2), 
##'                        M = 100)
##'                        
##' oldpar <- par(mfrow = c(2, 2))
##' plot(survfit, estimator = "both", include.y = TRUE)
##' par(oldpar)
##' }
##' @export
##' 
plot.survfitjmcs <- function (x, estimator = c("both", "mean", "median"), 
                            conf.int = TRUE, include.y = FALSE, ...) {
  
  if (!inherits(x, "survfitjmcs"))
    stop("Use only with 'survfitjmcs' xs.\n")
  
  if (!x$CompetingRisk) {
    which = 1:nrow(x$Last.time)
    ylim <- c(0, 1)
    
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
        probmedian <- c(1, 1, x$Pred[[i]][, 3])
        if (estimator == "both") {
          if (!include.y) {
            plot(times, probmean, xlab = "Time", ylab = expression(paste("Pr(", T[i] >= u, " | ", T[i] > s, 
                                                                         ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), 
                 main = paste("Subject", x$Last.time[i, 1], sep = " "), col = "red", type = "l", ylim = ylim)
            lines(times, probmedian, col = "green", type = "l")
            segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                     y1 = 1,
                     lwd = 1)
          } else {
            plot(x$y.obs[[i]][, 1], x$y.obs[[i]][, 2], xlim = c(0, max(x$Pred[[i]][, 1])), axes = TRUE, xlab = "Time", 
                 ylab = "Longitudinal outcome", type = "p", pch = 8)
            par(new = TRUE)    
            plot(times, probmean, xlab = "", ylab = "", 
                 main = paste("Subject", x$Last.time[i, 1], sep = " "), col = "red", type = "l", ylim = ylim, axes = FALSE)
            axis(side = 4, at = pretty(range(c(0, 1)))) 
            mtext(expression(paste("Pr(", T[i] >= u, " | ", T[i] > s, 
                                   ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), side=4, line=3)
            lines(times, probmedian, col = "green", type = "l")
            segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                     y1 = 1,
                     lwd = 1)
          }
        } else if (estimator == "median") {
          if (!include.y) {
            plot(times, probmedian, xlab = "Time", ylab = expression(paste("Pr(", T[i] >= u, " | ", T[i] > s, 
                                                                           ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), 
                 main = paste("Subject", x$Last.time[i, 1], sep = " "), col = "green", type = "l", ylim = ylim)
            abline(v = as.numeric(x$Last.time[i, 2]))
            segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                     y1 = 1,
                     lwd = 1)
          } else {
            plot(x$y.obs[[i]][, 1], x$y.obs[[i]][, 2], xlim = c(0, max(x$Pred[[i]][, 1])), axes = TRUE, xlab = "Time", 
                 ylab = "Longitudinal outcome", type = "p", pch = 8)
            par(new = TRUE)    
            plot(times, probmedian, xlab = "", ylab = "", 
                 main = paste("Subject", x$Last.time[i, 1], sep = " "), col = "green", type = "l", ylim = ylim, axes = FALSE)
            axis(side = 4, at = pretty(range(c(0, 1)))) 
            mtext(expression(paste("Pr(", T[i] >= u, " | ", T[i] > s, 
                                   ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), side=4, line=3)
            segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                     y1 = 1,
                     lwd = 1)
          }
        } else if (estimator == "mean") {
          if (!include.y) {
            plot(times, probmean, xlab = "Time", ylab = expression(paste("Pr(", T[i] >= u, " | ", T[i] > s, 
                                                                         ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), 
                 main = paste("Subject", x$Last.time[i, 1], sep = " "), col = "red", type = "l", ylim = ylim)
            abline(v = as.numeric(x$Last.time[i, 2]))
            segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                     y1 = 1,
                     lwd = 1)
          } else {
            plot(x$y.obs[[i]][, 1], x$y.obs[[i]][, 2], xlim = c(0, max(x$Pred[[i]][, 1])), axes = TRUE, xlab = "Time", 
                 ylab = "Longitudinal outcome", type = "p", pch = 8)
            par(new = TRUE)    
            plot(times, probmean, xlab = "", ylab = "", 
                 main = paste("Subject", x$Last.time[i, 1], sep = " "), col = "red", type = "l", ylim = ylim, axes = FALSE)
            axis(side = 4, at = pretty(range(c(0, 1)))) 
            mtext(expression(paste("Pr(", T[i] >= u, " | ", T[i] > s, 
                                   ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), side=4, line=3)
            segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                     y1 = 1,
                     lwd = 1)
          }
        } else {
          stop("Please choose one of the following estimators: mean, median, or both.")
        }
        
        if (conf.int) {
          Lower <- c(1, 1, x$Pred[[i]][, 4])
          Upper <- c(1, 1, x$Pred[[i]][, 5])
          lines(times, Lower, col = "black", type = "l", lty = 2)
          lines(times, Upper, col = "black", type = "l", lty = 2)
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
        if (show[(i-1)*2+j] && !return) {
          times <- c(0, as.numeric(x$Last.time[i, 2]), x$Pred[[i]][[j]][, 1])
          probmean <- c(0, 0, x$Pred[[i]][[j]][, 2])
          probmedian <- c(0, 0, x$Pred[[i]][[j]][, 3])
          if (estimator == "both") {
            if (!include.y) {
              plot(times, probmean, xlab = "Time", ylab = expression(paste("Pr(", T[i] <= u, ",", D[i] == k, " | ", T[i] > s, 
                                                                           ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), 
                   main = paste("Subject", x$Last.time[i, 1], "k=", j, sep = " "), col = "red", type = "l", ylim = ylim)
              lines(times, probmedian, col = "green", type = "l")
              segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                       y1 = 1,
                       lwd = 1)
            } else {
              plot(x$y.obs[[i]][, 1], x$y.obs[[i]][, 2], xlim = c(0, max(x$Pred[[i]][[j]][, 1])), axes = TRUE, xlab = "Time", 
                   ylab = "Longitudinal outcome", type = "p", pch = 8)
              par(new = TRUE)    
              plot(times, probmean, xlab = "", ylab = "", 
                   main = paste("Subject", x$Last.time[i, 1], "k=", j, sep = " "), col = "red", type = "l", ylim = ylim, axes = FALSE)
              axis(side = 4, at = pretty(range(c(0, 1)))) 
              mtext(expression(paste("Pr(", T[i] <= u, ", ", D[i] == k, " | ", T[i] > s, 
                                     ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), side=4, line=3)
              lines(times, probmedian, col = "green", type = "l")
              segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                       y1 = 1,
                       lwd = 1)
            }
          } else if (estimator == "median") {
            if (!include.y) {
              plot(times, probmedian, xlab = "Time", ylab = expression(paste("Pr(", T[i] <= u, ",", D[i] == k, " | ", T[i] > s, 
                                                                             ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), 
                   main = paste("Subject", x$Last.time[i, 1], "k=", j, sep = " "), col = "green", type = "l", ylim = ylim)
              abline(v = as.numeric(x$Last.time[i, 2]))
              segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                       y1 = 1,
                       lwd = 1)
            } else {
              plot(x$y.obs[[i]][, 1], x$y.obs[[i]][, 2], xlim = c(0, max(x$Pred[[i]][[j]][, 1])), axes = TRUE, xlab = "Time", 
                   ylab = "Longitudinal outcome", type = "p", pch = 8)
              par(new = TRUE)    
              plot(times, probmedian, xlab = "", ylab = "", 
                   main = paste("Subject", x$Last.time[i, 1], "k=", j, sep = " "), col = "green", type = "l", ylim = ylim, axes = FALSE)
              axis(side = 4, at = pretty(range(c(0, 1)))) 
              mtext(expression(paste("Pr(", T[i] <= u, ",", D[i] == k, " | ", T[i] > s, 
                                     ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), side=4, line=3)
              segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                       y1 = 1,
                       lwd = 1)
            }
          } else if (estimator == "mean") {
            if (!include.y) {
              plot(times, probmean, xlab = "Time", ylab = expression(paste("Pr(", T[i] <= u, ",", D[i] == k, " | ", T[i] > s, 
                                                                           ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), 
                   main = paste("Subject", x$Last.time[i, 1], "k=", j, sep = " "), col = "red", type = "l", ylim = ylim)
              abline(v = as.numeric(x$Last.time[i, 2]))
              segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                       y1 = 1,
                       lwd = 1)
            } else {
              plot(x$y.obs[[i]][, 1], x$y.obs[[i]][, 2], xlim = c(0, max(x$Pred[[i]][[j]][, 1])), axes = TRUE, xlab = "Time", 
                   ylab = "Longitudinal outcome", type = "p", pch = 8)
              par(new = TRUE)    
              plot(times, probmean, xlab = "", ylab = "", 
                   main = paste("Subject", x$Last.time[i, 1], "k=", j, sep = " "), col = "red", type = "l", ylim = ylim, axes = FALSE)
              axis(side = 4, at = pretty(range(c(0, 1)))) 
              mtext(expression(paste("Pr(", T[i] <= u, ", ", D[i] == k, " | ", T[i] > s, 
                                     ", ", y[i]^(s), ", ",  Psi,")", sep = " ")), side=4, line=3)
              segments(x0 = as.numeric(x$Last.time[i, 2]), x1 = as.numeric(x$Last.time[i, 2]), y0 = -1,
                       y1 = 1,
                       lwd = 1)
            }
          } else {
            stop("Please choose one of the following estimators: mean, median, or both.")
          }
          
          if (conf.int) {
            Lower <- c(0, 0, x$Pred[[i]][[j]][, 4])
            Upper <- c(0, 0, x$Pred[[i]][[j]][, 5])
            lines(times, Lower, col = "black", type = "l", lty = 2)
            lines(times, Upper, col = "black", type = "l", lty = 2)
          }
        }
      }
    }
    invisible()
  }
}