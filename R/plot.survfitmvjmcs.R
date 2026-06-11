##' @title Plot predictions from a multivariate joint model
##' @name plot
##' @aliases plot.survfitmvjmcs
##'
##' @description
##' Plot subject-specific conditional survival probabilities or cumulative
##' incidence probabilities from a \code{survfitmvjmcs} object. Longitudinal
##' observations for multiple biomarkers are stacked vertically, with the
##' predicted survival or cumulative incidence curve overlaid on each panel.
##'
##' @param x An object of class \code{survfitmvjmcs}.
##' @param subject Optional subject or subjects to plot. Can be a subject ID,
##' a numeric subject index, or a vector of IDs or indices. If \code{NULL}, all
##' subjects in \code{x$Last.time} are plotted. Default is \code{NULL}.
##' @param risk Failure type to plot when \code{x$CompetingRisk = TRUE}.
##' Default is \code{1}.
##' @param include.y Logical; retained for consistency with other plotting
##' methods. The current method always displays longitudinal biomarker panels.
##' Default is \code{TRUE}.
##' @param surv.panel A character to specify where the predicted
##' survival or cumulative incidence curve is displayed. Options are
##' \code{"first"} to overlay the curve on the first biomarker panel,
##' \code{"last"} to overlay it on the last biomarker panel, \code{"all"} to
##' overlay it on all biomarker panels. Default is \code{"first"}.
##' @param xlab X axis label.
##' @param ylab Y axis label.
##' @param xlim X axis support.
##' @param ylim.long Y axis support for the longitudinal outcome.
##' @param ylim.surv Y axis support for the event / survival probability.
##' @param ... Additional graphical arguments passed to the biomarker point plot.
##'
##' @return Invisibly returns \code{x}. The function is called for its plotting
##' side effect.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{survfitmvjmcs}}
##' @export

plot.survfitmvjmcs <- function(
    x,
    subject = NULL,
    risk = 1,
    include.y = TRUE,
    surv.panel = "first",
    xlab = "Time", 
    ylab = NULL, 
    xlim = NULL, 
    ylim.long = NULL, 
    ylim.surv = c(0, 1),
    ...
) {
  
  if (!inherits(x, "survfitmvjmcs")) {
    stop("Use only with 'survfitmvjmcs' objects.\n")
  }
  
  # Reset plotting overlay state in case a previous failed plot left par(new = TRUE).
  # This helps avoid warnings like:
  #   "calling par(new=TRUE) with no plot"
  par(new = FALSE)
  
  xlab <- "Time"
  xlim <- NULL
  ylim.long <- NULL
  ylim.surv <- c(0, 1)
  
  col.surv <- "red"
  lwd.surv <- 2.5
  
  col.points <- "black"
  pch <- 8
  cex.points <- 1.2
  
  lwd.landmark <- 1.2
  lwd.box <- 1.3
  lwd.axis <- 1.1
  
  cex.axis <- 0.9
  cex.lab <- 0.95
  cex.main <- 1.15
  
  subject.ids <- x$Last.time[, 1]
  subject.ids.char <- as.character(subject.ids)
  
  # If subject is NULL, plot every subject in x$Last.time.
  if (is.null(subject)) {
    
    subject.indices <- seq_along(subject.ids)
    
  } else {
    
    # First try matching by actual subject ID.
    subject.match <- match(as.character(subject), subject.ids.char)
    
    if (all(!is.na(subject.match))) {
      
      subject.indices <- subject.match
      
    } else if (
      is.numeric(subject) &&
      all(subject %in% seq_along(subject.ids))
    ) {
      
      # If matching by subject ID fails, allow numeric indexing.
      subject.indices <- as.integer(subject)
      
    } else {
      
      stop(
        "Requested subject(s) not found. Use subject IDs from x$Last.time[, 1], ",
        "or valid numeric indices."
      )
    }
  }
  
  
  numBio <- length(x$y.obs)
  
  if (numBio < 1) {
    stop("No longitudinal biomarker data found in x$y.obs.")
  }
  
  
  # Decide which stacked panel gets the red survival/CIF curve
  if (identical(surv.panel, "all")) {
    
    surv.panels <- seq_len(numBio)
    
  } else if (identical(surv.panel, "last")) {
    
    surv.panels <- numBio
    
  } else if (identical(surv.panel, "first")) {
    
    surv.panels <- 1
    
  } else if (is.numeric(surv.panel) && surv.panel %in% seq_len(numBio)) {
    
    surv.panels <- as.integer(surv.panel)
    
  } else {
    
    stop("surv.panel must be 'all', 'last', 'first', or a valid biomarker number.")
  }
  
  
  # Recover real biomarker names
  biomarker.names <- names(x$y.obs)
  
  if (
    is.null(biomarker.names) ||
    length(biomarker.names) != numBio ||
    any(is.na(biomarker.names)) ||
    any(biomarker.names == "")
  ) {
    
    first.i <- subject.indices[1]
    
    biomarker.names <- sapply(seq_len(numBio), function(g) {
      
      ydat.g <- x$y.obs[[g]][[first.i]]
      
      if (!is.null(ydat.g) && ncol(ydat.g) >= 2) {
        
        possible.name <- colnames(ydat.g)[2]
        
        if (
          !is.null(possible.name) &&
          !is.na(possible.name) &&
          possible.name != ""
        ) {
          return(possible.name)
        }
      }
      
      paste("Biomarker", g)
    })
  }
  
  ylab <- biomarker.names
  
  
  # Plot dimensions / graphical parameters
  oldpar <- par(c("mfrow", "mar", "oma", "mgp", "tcl", "ask"))
  
  on.exit({
    par(new = FALSE)
    do.call(par, oldpar)
  }, add = TRUE)
  
  
  # If plotting multiple subjects interactively, show one subject at a time.
  if (length(subject.indices) > 1 && interactive()) {
    par(ask = TRUE)
  }
  
  
  # Looping over subjects
  for (i in subject.indices) {
    
    subject.id <- x$Last.time[i, 1]
    last.time <- as.numeric(x$Last.time[i, 2])
    
    # Extract prediction curve for subject.
    pred.times.full <- x$Pred[[i]][, 1]
    
    # Keep prediction times only at or after the landmark time.
    #
    # This removes the artificial flat red line before the vertical landmark.
    # In the older version, the curve used:
    #   times <- c(0, last.time, pred.times)
    # which forced a red line from time 0 to last.time.
    #
    # Now the curve starts at last.time instead.
    keep.pred <- pred.times.full >= last.time
    pred.times <- pred.times.full[keep.pred]
    
    
    # Single failure type case
    
    if (!x$CompetingRisk) {
      
      # Single failure / survival probability case.
      pred.probs <- x$Pred[[i]][keep.pred, 2]
      
      # Start curve at landmark time.
      #
      # Survival probability is 1 at the landmark because the plotted quantity
      # is conditional on the subject surviving past the landmark time s.
      times <- c(last.time, pred.times)
      probmean <- c(1, pred.probs)
      
      prob.ylab <- expression(
        paste(
          "Pr(", T[i] >= u, " | ", T[i] > s,
          ", ", y[i]^(s), ", ", Psi, ")",
          sep = " "
        )
      )
      
      main.title <- paste("Subject", subject.id, "- Single failure type")
      
    } else {
      
      # Competing risk / cumulative incidence case.
      n.risks <- ncol(x$Pred[[i]]) - 1
      
      if (!risk %in% seq_len(n.risks)) {
        stop("Requested risk is not available in x$Pred[[subject]].")
      }
      
      pred.probs <- x$Pred[[i]][keep.pred, risk + 1]
      
      # Start curve at landmark time.
      #
      # CIF is 0 at the landmark because the plotted quantity is conditional
      # on no event occurring before landmark time s.
      times <- c(last.time, pred.times)
      probmean <- c(0, pred.probs)
      
      prob.ylab <- bquote(
        Pr(T[i] <= u, D[i] == .(risk) ~ "|" ~
             T[i] > s, y[i]^{(s)}, Psi)
      )
      
      main.title <- paste(
        "Subject",
        subject.id,
        "- Competing risks: risk",
        risk
      )
    }
    
    
    # x-axis specifications
    if (is.null(xlim)) {
      
      all.long.times <- c()
      
      for (g in seq_len(numBio)) {
        if (!is.null(x$y.obs[[g]][[i]])) {
          all.long.times <- c(all.long.times, x$y.obs[[g]][[i]][, 1])
        }
      }
      
      xlim.i <- c(
        0,
        max(c(pred.times.full, all.long.times), na.rm = TRUE)
      )
      
    } else {
      
      xlim.i <- xlim
    }
    
    
    # Stacked layout for each subject
    par(
      mfrow = c(numBio, 1),
      oma = c(2.6, 0, 2.0, 3.2),
      mgp = c(1.8, 0.5, 0),
      tcl = -0.25
    )
    
    
    # Plot one biomarker per row
    for (g in seq_len(numBio)) {
      
      ydat <- x$y.obs[[g]][[i]]
      
      if (is.null(ydat) || nrow(ydat) == 0) {
        warning(
          paste(
            "No longitudinal observations found for biomarker",
            g,
            "and subject",
            subject.id
          )
        )
        next
      }
      
      ytime <- ydat[, 1]
      yvalue <- ydat[, 2]
      
      
      # Biomarker y-axis
      if (is.null(ylim.long)) {
        
        ylim.g <- range(yvalue, na.rm = TRUE)
        
        if (diff(ylim.g) == 0) {
          ylim.g <- ylim.g + c(-0.5, 0.5)
        }
        
      } else if (is.list(ylim.long)) {
        
        ylim.g <- ylim.long[[g]]
        
      } else {
        
        ylim.g <- ylim.long
      }
      
      
      # PANEL MARGINS
      if (numBio == 1) {
        
        par(mar = c(3.0, 3.8, 1.5, 3.2))
        
      } else if (g == 1) {
        
        # Top panel: small bottom margin
        par(mar = c(0.3, 3.8, 0.9, 3.2))
        
      } else if (g == numBio) {
        
        # Bottom panel: small top margin
        par(mar = c(2.4, 3.8, 0.3, 3.2))
        
      } else {
        
        # Middle panels: minimal vertical margin
        par(mar = c(0.3, 3.8, 0.3, 3.2))
      }
      
      
      # Main biomarker plot
      plot(
        ytime,
        yvalue,
        xlim = xlim.i,
        ylim = ylim.g,
        xlab = "",
        ylab = "",
        type = "p",
        pch = pch,
        cex = cex.points,
        col = col.points,
        axes = FALSE,
        ...
      )
      
      
      # Left biomarker axis
      axis(
        side = 2,
        lwd = lwd.axis,
        lwd.ticks = lwd.axis,
        cex.axis = cex.axis
      )
      
      
      # Left biomarker label
      mtext(
        ylab[g],
        side = 2,
        line = 2.4,
        cex = cex.lab
      )
      
      
      # Bottom time axis only on final biomarker panel
      if (g == numBio) {
        axis(
          side = 1,
          lwd = lwd.axis,
          lwd.ticks = lwd.axis,
          cex.axis = cex.axis
        )
      }
      
      
      # Border before overlay
      box(lwd = lwd.box)
      
      
      # Landmark line
      segments(
        x0 = last.time,
        x1 = last.time,
        y0 = ylim.g[1],
        y1 = ylim.g[2],
        lwd = lwd.landmark,
        col = "black"
      )
      
      
      # Survival/CIF overlay
      if (g %in% surv.panels) {
        
        # Plotting survival/CIF curve.
        par(new = TRUE)
        
        plot(
          times,
          probmean,
          xlim = xlim.i,
          ylim = ylim.surv,
          xlab = "",
          ylab = "",
          type = "l",
          col = col.surv,
          lwd = lwd.surv,
          axes = FALSE
        )
        
        # Right probability axis for this panel.
        axis(
          side = 4,
          at = pretty(ylim.surv),
          las = 2,
          lwd = lwd.axis,
          lwd.ticks = lwd.axis,
          cex.axis = cex.axis
        )
        
        # Redraw border after overlay.
        box(lwd = lwd.box)
      }
    }
    
    
    # Labels
    mtext(
      main.title,
      side = 3,
      outer = TRUE,
      line = 0.5,
      font = 2,
      cex = cex.main
    )
    
    mtext(
      xlab,
      side = 1,
      outer = TRUE,
      line = 1.1,
      cex = cex.lab
    )
    
    # Centered right-side probability label across all biomarker panels.
    mtext(
      prob.ylab,
      side = 4,
      outer = TRUE,
      line = 2.0,
      cex = cex.lab
    )
  }
  
  invisible(x)
}
