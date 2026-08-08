##' @title Plot conditional probabilities for new subjects
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
##' @param xlab X axis label.
##' @param ylab Y axis label.
##' @param xlim X axis support.
##' @param ylim.surv Y axis support for the event / survival probability.
##' @param ... Additional graphical arguments passed to the biomarker point plot.
##'
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{survfitmvjmcs}}
##' @export

plot.survfitmvjmcs <- function(x,
                               subject = NULL,
                               risk = 1,
                               include.y = TRUE,
                               xlab = NULL, ylab = NULL, 
                               xlim = NULL,
                               ylim.surv = NULL, 
                               ...) {
  
  if (!inherits(x, "survfitmvjmcs")) {
    stop("Use only with 'survfitmvjmcs' objects.\n")
  }
  
  # Reset plotting overlay state in case a previous failed plot left par(new = TRUE)
  # Helps avoid warnings like:
  #   "calling par(new=TRUE) with no plot"
  par(new = FALSE)
  
  if (is.null(xlab)) xlab <- "Time"
  if (is.null(ylim.surv)) ylim.surv <- c(0, 1)
  
  col.surv <- "red"
  lwd.surv <- 2.5
  
  col.points <- "black"
  pch.points <- 16
  cex.points <- 0.75
  
  lwd.landmark <- 1.2
  lwd.box <- 1.3
  lwd.axis <- 1.1
  
  cex.axis <- 0.9
  cex.lab <- 0.95
  cex.main <- 1.15
  
  subject.ids <- x$Last.time[, 1]
  subject.ids.char <- as.character(subject.ids)
  
  # If subject is NULL, plot every subject in x$Last.time
  if (is.null(subject)) {
    
    subject.indices <- seq_along(subject.ids)
    
  } else {
    
    # First try matching by actual subject ID
    subject.match <- match(as.character(subject), subject.ids.char)
    
    if (all(!is.na(subject.match))) {
      
      subject.indices <- subject.match
      
    } else if (
      is.numeric(subject) &&
      all(subject %in% seq_along(subject.ids))
    ) {
      
      # If matching by subject ID fails, allow numeric indexing
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
  
  if (!x$CompetingRisk) {
    stop("This panel plot is intended for competing-risk survfitmvjmcs objects.")
  }
  
  if (length(subject.indices) < 2) {
    stop("Need at least two subjects, e.g., subject = c(2432, 4157).")
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
  
  if (is.null(ylab)) ylab <- biomarker.names
  
  # Plot dimensions / graphical parameters
  oldpar <- par(c("mfrow", "mar", "oma", "mgp", "tcl", "ask"))
  
  on.exit({
    par(new = FALSE)
    do.call(par, oldpar)
  }, add = TRUE)
  
  n.subjects.panel <- length(subject.indices)
  n.rows.panel <- numBio + 1
  
  # Use one page for all selected subjects
  # Rows = biomarkers + CIF row
  # Columns = selected subjects
  par(
    mfrow = c(n.rows.panel, n.subjects.panel),
    oma = c(3.0, 3.8, 3.0, 1.5),
    mgp = c(1.8, 0.5, 0),
    tcl = -0.25,
    ask = FALSE
  )
  
  # Same x-axis range across all panels
  all.panel.times <- c()
  
  for (ii in subject.indices) {
    
    all.panel.times <- c(all.panel.times, x$Pred[[ii]][, 1])
    
    for (g in seq_len(numBio)) {
      if (!is.null(x$y.obs[[g]][[ii]])) {
        all.panel.times <- c(all.panel.times, x$y.obs[[g]][[ii]][, 1])
      }
    }
  }
  
  xlim <- c(0, max(all.panel.times, na.rm = TRUE))
  
  # Same y-axis range for each biomarker across all selected subjects
  # This makes side-by-side subject comparisons easier
  ylim.panel.long <- vector("list", numBio)
  
  for (g in seq_len(numBio)) {
    
    all.y.g <- c()
    
    for (ii in subject.indices) {
      
      ydat.g <- x$y.obs[[g]][[ii]]
      
      if (!is.null(ydat.g) && nrow(ydat.g) > 0) {
        all.y.g <- c(all.y.g, ydat.g[, 2])
      }
    }
    
    ylim.g <- range(all.y.g, na.rm = TRUE)
    
    if (diff(ylim.g) == 0) {
      ylim.g <- ylim.g + c(-0.5, 0.5)
    }
    
    ylim.panel.long[[g]] <- ylim.g
  }
  
  # Plot biomarker rows.
  for (g in seq_len(numBio)) {
    
    for (col.i in seq_along(subject.indices)) {
      
      ii <- subject.indices[col.i]
      
      subject.id.panel <- x$Last.time[ii, 1]
      last.time.panel <- as.numeric(x$Last.time[ii, 2])
      
      ydat <- x$y.obs[[g]][[ii]]
      
      # PANEL MARGINS
      if (g == 1) {
        
        # Top biomarker row.
        par(mar = c(0.3, 3.2, 1.4, 0.6))
        
      } else {
        
        # Middle biomarker rows.
        par(mar = c(0.3, 3.2, 0.3, 0.6))
      }
      
      if (is.null(ydat) || nrow(ydat) == 0) {
        
        plot(
          NA,
          NA,
          xlim = xlim,
          ylim = ylim.panel.long[[g]],
          xlab = "",
          ylab = "",
          axes = FALSE
        )
        
        warning(
          paste(
            "No longitudinal observations found for biomarker",
            g,
            "and ID",
            subject.id.panel
          )
        )
        
      } else {
        
        ytime <- ydat[, 1]
        yvalue <- ydat[, 2]
        
        plot(
          ytime,
          yvalue,
          xlim = xlim,
          ylim = ylim.panel.long[[g]],
          xlab = "",
          ylab = "",
          type = "b",
          pch = pch.points,
          cex = cex.points,
          lwd = 1.1,
          col = col.points,
          axes = FALSE,
          ...
        )
      }
      
      # Left y-axis only for the first column
      if (col.i == 1) {
        axis(
          side = 2,
          lwd = lwd.axis,
          lwd.ticks = lwd.axis,
          cex.axis = cex.axis
        )
        
        mtext(
          ylab[g],
          side = 2,
          line = 2.3,
          cex = cex.lab
        )
      }
      
      # Top title only on first biomarker row.
      if (g == 1) {
        title(
          main = paste("ID", subject.id.panel),
          cex.main = 0.95,
          line = 0.4
        )
      }
      
      box(lwd = lwd.box)
      
      # Landmark line.
      segments(
        x0 = last.time.panel,
        x1 = last.time.panel,
        y0 = ylim.panel.long[[g]][1],
        y1 = ylim.panel.long[[g]][2],
        lwd = lwd.landmark,
        col = "black"
      )
    }
  }
  
  # Plot bottom CIF row
  for (col.i in seq_along(subject.indices)) {
    
    ii <- subject.indices[col.i]
    
    last.time.panel <- as.numeric(x$Last.time[ii, 2])
    pred.times.full.panel <- x$Pred[[ii]][, 1]
    
    keep.pred.panel <- pred.times.full.panel >= last.time.panel
    pred.times.panel <- pred.times.full.panel[keep.pred.panel]
    
    n.risks.panel <- ncol(x$Pred[[ii]]) - 1
    
    if (!risk %in% seq_len(n.risks.panel)) {
      stop("Requested risk is not available in x$Pred[[subject]].")
    }
    
    pred.probs.panel <- x$Pred[[ii]][keep.pred.panel, risk + 1]
    
    # Start curve at landmark time.
    #
    # CIF is 0 at the landmark because the plotted quantity is conditional
    # on no event occurring before landmark time s
    times.panel <- c(last.time.panel, pred.times.panel)
    probmean.panel <- c(0, pred.probs.panel)
    
    par(mar = c(2.4, 3.2, 0.3, 0.6))
    
    plot(
      times.panel,
      probmean.panel,
      xlim = xlim,
      ylim = ylim.surv,
      xlab = "",
      ylab = "",
      type = "l",
      col = col.surv,
      lwd = lwd.surv,
      axes = FALSE
    )
    
    # Left y-axis only for first column.
    if (col.i == 1) {
      axis(
        side = 2,
        at = pretty(ylim.surv),
        lwd = lwd.axis,
        lwd.ticks = lwd.axis,
        cex.axis = cex.axis
      )
      
      mtext(
        "CIF",
        side = 2,
        line = 2.3,
        cex = cex.lab
      )
    }
    
    # Bottom x-axis for the CIF row
    axis(
      side = 1,
      lwd = lwd.axis,
      lwd.ticks = lwd.axis,
      cex.axis = cex.axis
    )
    
    box(lwd = lwd.box)
    
    # Landmark line
    segments(
      x0 = last.time.panel,
      x1 = last.time.panel,
      y0 = ylim.surv[1],
      y1 = ylim.surv[2],
      lwd = lwd.landmark,
      col = "black"
    )
  }
  
  # Labels.
  mtext(paste("Competing risks: risk", risk),
        side = 3,
        outer = TRUE,
        line = 1.0,
        font = 2,
        cex = cex.main)
  
  mtext(xlab,
        side = 1,
        outer = TRUE,
        line = 1.1,
        cex = cex.lab)
  
  invisible(x)
}