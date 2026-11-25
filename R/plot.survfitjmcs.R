##' @title Plot dynamic predictions from joint models
##' @name plot.survfit
##' @description Plot subject-specific dynamic predictions, optionally
##'   together with observed longitudinal responses.
##'
##' @param x An object of class \code{survfitjmcs},
##'   \code{survfitJMMLSM}, or \code{survfitmvjmcs}.
##' @param include.y Logical. For \code{survfitjmcs} and
##'   \code{survfitJMMLSM}, include observed longitudinal responses
##'   when \code{TRUE}; the default is \code{FALSE}.
##'   For \code{survfitmvjmcs}, the default is \code{TRUE}, and
##'   longitudinal biomarker panels are always displayed regardless
##'   of this argument.
##' @param xlab Label for the time axis. If \code{NULL},
##'   \code{"Time"} is used.
##' @param ylab Label for the response or probability axis, as
##'   applicable. For \code{survfitmvjmcs}, a character vector with
##'   one label per biomarker. If \code{NULL}, default labels are used.
##' @param xlim Numeric vector of length two specifying the time-axis
##'   limits, or \code{NULL} for automatic limits.
##' @param ylim.long Numeric vector of length two specifying the
##'   longitudinal response-axis limits. Applies to
##'   \code{survfitjmcs} and \code{survfitJMMLSM} when
##'   \code{include.y = TRUE}.
##' @param ylim.surv Numeric vector of length two specifying the
##'   probability-axis limits. If \code{NULL}, method-specific
##'   defaults are used.
##' @param subject For \code{survfitmvjmcs} only, subject IDs or
##'   numeric subject indices. IDs are matched first; numeric indices
##'   are used if ID matching fails. If \code{NULL}, all subjects in
##'   \code{x$Last.time} are selected. At least two subjects must
##'   be selected.
##' @param risk For \code{survfitmvjmcs} only, the event type whose
##'   cumulative incidence is plotted. Default is \code{1}.
##' @param ... For \code{survfitmvjmcs}, additional graphical
##'   arguments passed to the biomarker point plots.
##'   Currently unused by the other two methods.
##'
##' @details
##' For \code{survfitjmcs} and \code{survfitJMMLSM}, single-failure
##' models produce conditional survival curves, and competing-risks
##' models produce conditional cumulative incidence curves.
##'
##' The current \code{survfitmvjmcs} method supports competing-risks
##' predictions only. Subjects are arranged in columns, biomarker
##' trajectories in the upper rows, and cumulative incidence curves
##' for the selected event type in the bottom row.
##'
##' @return Called for its graphical side effects.
##'   The \code{survfitjmcs} and \code{survfitJMMLSM} methods return
##'   \code{NULL} invisibly. The \code{survfitmvjmcs} method returns
##'   \code{x} invisibly.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{survfitJM}}
##' @export

plot.survfitjmcs <- function(
    x,
    include.y = FALSE,
    xlab = NULL,
    ylab = NULL,
    xlim = NULL,
    ylim.long = NULL,
    ylim.surv = NULL,
    ...
) {
  
  if (!inherits(x, "survfitjmcs"))
    stop("Use only with 'survfitjmcs' objects.\n")
  
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
        
        times <- c(as.numeric(x$Last.time[i, 2]), x$Pred[[i]][, 1])
        probmean <- c(1, x$Pred[[i]][, 2])
        
        if (is.null(xlim)) xlim <- c(0, max(x$Pred[[i]][, 1]))
        
        if (!include.y) {
          
          if (is.null(ylab)) {
            ylab <- expression(
              paste(
                "Pr(", T[i] >= u, " | ", T[i] > s,
                ", ", y[i]^(s), ", ", Psi, ")",
                sep = " "
              )
            )
          }
          
          plot(
            times,
            probmean,
            xlab = xlab,
            ylab = ylab,
            xlim = xlim,
            main = paste("ID", x$Last.time[i, 1], sep = " "),
            col = "red",
            type = "l",
            ylim = ylim.surv
          )
          
          segments(
            x0 = as.numeric(x$Last.time[i, 2]),
            x1 = as.numeric(x$Last.time[i, 2]),
            y0 = -1,
            y1 = 1,
            lwd = 1
          )
          
        } else {
          
          if (is.null(ylab)) ylab = "Longitudinal outcome"
          
          plot(
            x$y.obs[[i]][, 1],
            x$y.obs[[i]][, 2],
            xlim = xlim,
            axes = TRUE,
            xlab = xlab,
            ylab = "",
            type = "b",
            pch = 8,
            ylim = ylim.long
          )
          
          title(ylab = ylab, line = 2.5)
          
          par(new = TRUE)
          
          plot(
            times,
            probmean,
            xlab = "",
            ylab = "",
            main = paste("ID", x$Last.time[i, 1], sep = " "),
            xlim = xlim,
            col = "red",
            type = "l",
            ylim = ylim.surv,
            axes = FALSE
          )
          
          axis(side = 4, at = pretty(range(ylim.surv)), line = 0)
          
          mtext(
            expression(
              paste(
                "Pr(", T[i] >= u, " | ", T[i] > s,
                ", ", y[i]^(s), ", ", Psi, ")",
                sep = " "
              )
            ),
            side = 4,
            line = 2.5
          )
          
          segments(
            x0 = as.numeric(x$Last.time[i, 2]),
            x1 = as.numeric(x$Last.time[i, 2]),
            y0 = -1,
            y1 = 1,
            lwd = 1
          )
        }
      }
    }
    
    invisible()
    
  } else {
    
    which = 1:nrow(x$Last.time)
    ylim <- c(0, 1)
    
    ask = (prod(par("mfcol")) < length(which))
    show <- rep(TRUE, 2 * length(which))
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
        
        if (show[(i - 1) * 2 + j] && !return) {
          
          times <- c(as.numeric(x$Last.time[i, 2]), x$Pred[[i]][, 1])
          probmean <- c(0, x$Pred[[i]][, j + 1])
          
          if (!include.y) {
            
            plot(
              times,
              probmean,
              xlab = xlab,
              ylab = bquote(
                Pr(T[i] <= u, D[i] == .(j) ~ "|" ~ T[i] > s, ~ y[i]^(s), ~ Psi)
              ),
              main = paste("ID", x$Last.time[i, 1], "- Competing risks: risk", j, sep = " "),
              col = "red",
              type = "l",
              ylim = ylim.surv,
              xlim = xlim
            )
            
            segments(
              x0 = as.numeric(x$Last.time[i, 2]),
              x1 = as.numeric(x$Last.time[i, 2]),
              y0 = -1,
              y1 = 1,
              lwd = 1
            )
            
          } else {
            
            if (is.null(ylab)) ylab = "Longitudinal outcome"
            
            plot(
              x$y.obs[[i]][, 1],
              x$y.obs[[i]][, 2],
              xlim = xlim,
              axes = TRUE,
              xlab = xlab,
              ylab = "",
              type = "b",
              pch = 8,
              ylim = ylim.long
            )
            
            title(ylab = ylab, line = 2.5)
            
            par(new = TRUE)
            
            plot(
              times,
              probmean,
              xlab = "",
              ylab = "",
              main = paste("ID", x$Last.time[i, 1], "- Competing risks: risk", j, sep = " "),
              col = "red",
              type = "l",
              ylim = ylim.surv,
              axes = FALSE,
              xlim = xlim
            )
            
            axis(side = 4, at = pretty(range(ylim.surv)), line = 0)
            
            mtext(
              bquote(
                Pr(T[i] <= u, D[i] == .(j) ~ "|" ~ T[i] > s, ~ y[i]^(s), ~ Psi)
              ),
              side = 4,
              line = 3.5
            )
            
            segments(
              x0 = as.numeric(x$Last.time[i, 2]),
              x1 = as.numeric(x$Last.time[i, 2]),
              y0 = -1,
              y1 = 1,
              lwd = 1
            )
          }
        }
      }
    }
    
    invisible()
  }
}