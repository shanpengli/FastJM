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
    ...
) {
  
  if (!inherits(x, "survfitmvjmcs")) {
    stop("Use only with 'survfitmvjmcs' objects.\n")
  }
  
  # Plot specifications
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
  
  cex.axis <- 1
  cex.lab <- 1
  cex.main <- 1.2
  
  
  # Grabbing all subjects
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
  
  # Number of biomarkers
  numBio <- length(x$y.obs)
  
  if (numBio < 1) {
    stop("No longitudinal biomarker data found in x$y.obs.")
  }
  
  # Try to recover real biomarker names.
  #
  # Priority:
  #   1. names(x$y.obs)
  #   2. second column name of x$y.obs[[g]][[i]]
  #   3. fallback: "Biomarker 1", "Biomarker 2", etc. if 
  #                 no official name
  
  # Identifying biomarker names
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
  
  # Plot dimensions
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # If plotting multiple subjects interactively, show one subject at a time.
  if (length(subject.indices) > 1 && interactive()) {
    par(ask = TRUE)
  }
  
  # Looping over subjects
  for (i in subject.indices) {
    
    subject.id <- x$Last.time[i, 1]
    last.time <- as.numeric(x$Last.time[i, 2])
    
    # Extract prediction curve for subject
    pred.times <- x$Pred[[i]][, 1]
    
    if (!x$CompetingRisk) {
      
      # Single failure / survival probability case
      times <- c(0, last.time, pred.times)
      probmean <- c(1, 1, x$Pred[[i]][, 2])
      
      prob.ylab <- expression(
        paste(
          "Pr(", T[i] >= u, " | ", T[i] > s,
          ", ", y[i]^(s), ", ", Psi, ")",
          sep = " "
        )
      )
      
      main.title <- paste("Subject", subject.id)
      
    } else {
      
      # Competing risk / cumulative incidence case
      n.risks <- ncol(x$Pred[[i]]) - 1
      
      if (!risk %in% seq_len(n.risks)) {
        stop("Requested risk is not available in x$Pred[[subject]].")
      }
      
      times <- c(0, last.time, pred.times)
      probmean <- c(0, 0, x$Pred[[i]][, risk + 1])
      
      prob.ylab <- expression(
        paste(
          "Pr(", T[i] <= u, ",", D[i] == k, " | ",
          T[i] > s, ", ", y[i]^(s), ", ", Psi, ")",
          sep = " "
        )
      )
      
      main.title <- paste("Subject", subject.id, "k =", risk)
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
        max(c(pred.times, all.long.times), na.rm = TRUE)
      )
      
    } else {
      
      xlim.i <- xlim
    }
    
    # Stacked layout for each subject
    par(
      mfrow = c(numBio, 1),
      oma = c(3.2, 0, 2.4, 4.2),
      mgp = c(2.2, 0.7, 0),
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
      # mar = c(bottom, left, top, right)
      # No + 0.1 because that creates a visible internal gap.
      
      if (numBio == 1) {
        
        par(mar = c(3.5, 4.5, 2, 4.5))
        
      } else if (g == 1) {
        
        # Top panel: no bottom margin
        par(mar = c(0, 4.5, 1.2, 4.5))
        
      } else if (g == numBio) {
        
        # Bottom panel: no top margin
        par(mar = c(2.6, 4.5, 0, 4.5))
        
      } else {
        
        # Middle panels: no vertical margin
        par(mar = c(0, 4.5, 0, 4.5))
      }
      
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
        line = 2.8,
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
      
      # Plotting survival/CIF curve
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
      
      # Right probability axis for this panel
      axis(
        side = 4,
        at = pretty(ylim.surv),
        las = 2,
        lwd = lwd.axis,
        lwd.ticks = lwd.axis,
        cex.axis = cex.axis
      )
      
      # Redraw border after overlay
      box(lwd = lwd.box)
    }
    
    
    # Labels
    mtext(
      main.title,
      side = 3,
      outer = TRUE,
      line = 0.6,
      font = 2,
      cex = cex.main
    )
    
    mtext(
      xlab,
      side = 1,
      outer = TRUE,
      line = 1.2,
      cex = cex.lab
    )
    
    # Centered right-side probability label across all biomarker panels
    mtext(
      prob.ylab,
      side = 4,
      outer = TRUE,
      line = 2.4,
      cex = cex.lab
    )
  }
  
  invisible(x)
}