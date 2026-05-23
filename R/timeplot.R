##' @title Diagnostic plots for the fitted joint model
##' @name timeplot
##' @aliases timeplot
##'
##' @description
##' Creates diagnostic plots for a fitted \code{jmcs} object from the
##' \pkg{FastJM} package, specifically longitudinal biomarker trajectory, 
##' log residual variance, and event curve. The function always produces:
##' \enumerate{
##'   \item a longitudinal biomarker trajectory plot, and
##'   \item a plot of log residual variance over follow-up time.
##' }
##'
##' It also produces one event plot based on \code{object$CompetingRisk}:
##' \itemize{
##'   \item if \code{object$CompetingRisk = FALSE}, a Kaplan-Meier plot is shown;
##'   \item if \code{object$CompetingRisk = TRUE}, a cumulative incidence plot is shown.
##' }
##'
##' The event time and status variables are inferred from
##' \code{object$SurvivalSubmodel} when \code{event_time_col} and
##' \code{event_status_col} are not supplied.
##'
##' @param object A fitted object of class \code{jmcs}.
##' @param biomarker Unquoted name of the longitudinal biomarker variable in
##'   \code{object$ydata}.
##' @param id_col Unquoted name of the subject identifier variable in
##'   \code{object$ydata}.
##' @param time_col Unquoted name of the longitudinal follow-up time variable in
##'   \code{object$ydata}. This variable is used on the x-axis for the
##'   biomarker trajectory and log residual variance plots.
##' @param visit_col Optional unquoted name of a visit variable in
##'   \code{object$ydata}. If supplied and found, visit-level summaries are
##'   used for the longitudinal diagnostic plots. If omitted, or if the supplied
##'   column is not found, derived time bins are used instead.
##' @param n.obs Numeric value specifying the number of subjects randomly selected
##'   for the longitudinal biomarker trajectory plot. Default is \code{100}.
##' @param seed Optional integer random seed used for reproducible subject
##'   sampling in the longitudinal biomarker trajectory plot. Default is
##'   \code{100}. If \code{NULL}, no seed is set.
##' @param window_days Numeric value giving the half-width of the time window used
##'   when selecting observations around each visit- or bin-level mean time.
##'   Default is \code{50}.
##' @param time_bin_width Numeric value giving the width of derived time bins when
##'   \code{visit_col} is not supplied or not found. Default is \code{180}. For
##'   small-scale time variables, such as time ranging from 0 to 5, use a smaller
##'   value such as \code{0.5} or \code{1}.
##' @param biomarker_y_lab Optional character string for the y-axis label of the
##'   biomarker trajectory plot. If \code{NULL}, the name of \code{biomarker}
##'   is used.
##' @param logvar_y_lab Character string for the y-axis label of the log residual
##'   variance plot. Default is \code{"Log Residual Variance"}.
##' @param event_x_lab Optional character string for the x-axis label of the event
##'   plot. If \code{NULL}, the name of the event-time variable is used.
##' @param event_y_lab Optional character string for the y-axis label of the event
##'   plot. If \code{NULL}, the default label from the Kaplan-Meier or cumulative
##'   incidence plotting function is used.
##' @param x_lab Optional character string for the x-axis label of the longitudinal
##'   plots. If \code{NULL}, the name of \code{time_col} is used.
##' @param x_break_by Numeric value controlling spacing of x-axis tick marks in
##'   the longitudinal plots. Default is \code{NULL}.
##' @param ylim_mean Optional numeric vector of length 2 giving y-axis limits for
##'   the biomarker trajectory plot.
##' @param ylim_logvar Optional numeric vector of length 2 giving y-axis limits for
##'   the log residual variance plot.
##' @param event_time_col Optional character string giving the event-time variable
##'   in \code{object$cdata}. If \code{NULL}, the function attempts to infer
##'   it as the first variable in \code{object$SurvivalSubmodel}, i.e. the
##'   first argument of \code{Surv(time, status)}.
##' @param event_status_col Optional character string giving the event-status
##'   variable in \code{object$cdata}. If \code{NULL}, the function attempts
##'   to infer it as the second variable in \code{object$SurvivalSubmodel},
##'   i.e. the second argument of \code{Surv(time, status)}.
##' @param fail_code Numeric or integer code in \code{event_status_col} denoting
##'   the event of interest. For single-failure models, this is the event code used
##'   for the Kaplan-Meier plot. For competing-risk models, this is the primary
##'   event code used in the cumulative incidence plot. If the model is
##'   single-failure and there is exactly one non-censoring status code,
##'   \code{fail_code} can be inferred.
##' @param cr_code Numeric or integer code in \code{event_status_col} denoting the
##'   competing event. Required for competing-risk models unless there is exactly
##'   one other non-censoring status code besides \code{fail_code}, in which case
##'   it can be inferred. Ignored for single-failure Kaplan-Meier plots.
##' @param censor_code Numeric or integer code in \code{event_status_col} denoting
##'   censoring. Default is \code{0}.
##' @param primary_event_label Character string used in event plot labels for the
##'   primary event. Default is \code{"Primary event"}.
##' @param competing_event_label Character string used in cumulative incidence plot
##'   labels for the competing event. Required for informative competing-risk
##'   labeling. Ignored for single-failure Kaplan-Meier plots. Default is
##'   \code{NULL}.
##' @param km_title Optional character string for the Kaplan-Meier plot title. If
##'   \code{NULL}, a default title is constructed from \code{primary_event_label}.
##' @param cif_title Optional character string for the cumulative incidence plot
##'   title. If \code{NULL}, a default title is constructed from
##'   \code{primary_event_label} and \code{competing_event_label}.
##' @param center_event_plot Logical; if \code{TRUE}, the three-plot layout places
##'   the two longitudinal diagnostic plots on the first row and centers the event
##'   plot on the second row. If \code{FALSE}, the three plots are stacked in one
##'   column. Default is \code{TRUE}.
##'
##' @details
##' The longitudinal trajectory plot displays individual biomarker trajectories for
##' a random subset of subjects, along with a summary mean curve. If
##' \code{visit_col} is supplied and found in \code{object$ydata}, summaries
##' are calculated by observed visit. Otherwise, the function creates derived
##' time bins using \code{floor(time_col / time_bin_width)}.
##'
##' The residual variance plot displays the log residual variance over visit- or
##' bin-level follow-up time. Subject-level fitted values are obtained using
##' \code{fitted(object, type = "Subject", process = "Longitudinal")}.
##'
##' The event plot is chosen automatically from the fitted model:
##' \itemize{
##'   \item For \code{CompetingRisk = FALSE}, the function creates a Kaplan-Meier
##'   plot using \code{fail_code} as the event code and \code{censor_code} as the
##'   censoring code.
##'   \item For \code{CompetingRisk = TRUE}, the function creates a cumulative
##'   incidence plot using \code{fail_code} as the primary event and
##'   \code{cr_code} as the competing event.
##' }
##'
##' In competing-risk settings, the model identifies event codes, but it does not
##' determine which event is scientifically "primary." The user should choose
##' \code{fail_code} based on the scientific question. For example, if
##' \code{status = 1} is heart failure and \code{status = 2} is death, then
##' \code{fail_code = 1} and \code{cr_code = 2} answer the question: "What is the
##' cumulative incidence of heart failure over time, accounting for death before
##' heart failure?"
##'
##' @return
##' Invisibly returns a list with components:
##' \describe{
##'   \item{\code{combined}}{The combined plot object.}
##'   \item{\code{p1}}{The biomarker trajectory plot.}
##'   \item{\code{p2}}{The log residual variance plot.}
##'   \item{\code{p_cif}}{The cumulative incidence plot if
##'     \code{object$CompetingRisk = TRUE}; otherwise \code{NULL}.}
##'   \item{\code{p_km}}{The Kaplan-Meier plot if
##'     \code{object$CompetingRisk = FALSE}; otherwise \code{NULL}.}
##'   \item{\code{cif_fit}}{The fitted cumulative incidence object if applicable;
##'     otherwise \code{NULL}.}
##'   \item{\code{km_fit}}{The fitted Kaplan-Meier object if applicable; otherwise
##'     \code{NULL}.}
##'   \item{\code{event_time_col}}{The event-time column used.}
##'   \item{\code{event_status_col}}{The event-status column used.}
##'   \item{\code{fail_code}}{The primary event code used.}
##'   \item{\code{cr_code}}{The competing event code used, or \code{NULL} for
##'     single-failure models.}
##'   \item{\code{censor_code}}{The censoring code used.}
##'   \item{\code{summary_data}}{A data frame of visit- or bin-level longitudinal
##'     summary statistics.}
##'   \item{\code{obs_per_subject_visit}}{A data frame giving the number of
##'     observations per subject and visit/bin after window-based filtering.}
##'   \item{\code{grouping_used}}{A character string indicating whether summaries
##'     were based on an observed visit variable or derived time bins.}
##' }
##'
##' @seealso
##' \code{\link[FastJM]{jmcs}}
##' 
##' @examples
##' \donttest{
##' library(FastJM)
##'
##' data(ydata)
##' data(cdata)
##'
##' # Single-failure example:
##' # Treat failure_type == 1 as the only event and all other observations as censored.
##' cdata_single <- cdata
##' cdata_single$event_single <- as.integer(cdata_single$failure_type == 1)
##'
##' fit_single <- jmcs(
##'   ydata = ydata,
##'   cdata = cdata_single,
##'   long.formula = response ~ time + gender + x1 + race,
##'   surv.formula = Surv(surv, event_single) ~ x1 + gender + x2 + race,
##'   random = ~ time | ID
##' )
##'
##' singleplot <- timeplot(
##'   object = fit_single,
##'   biomarker = response,
##'   id_col = ID,
##'   time_col = time,
##'   time_bin_width = 0.5,
##'   primary_event_label = "Event type 1"
##' )
##'
##' singleplot$combined
##' singleplot$p_km
##'
##' # Competing-risk example using FastJM simulated data.
##' fit_cr <- jmcs(
##'   ydata = ydata,
##'   cdata = cdata,
##'   long.formula = response ~ time + gender + x1 + race,
##'   surv.formula = Surv(surv, failure_type) ~ x1 + gender + x2 + race,
##'   random = ~ time | ID
##' )
##'
##' crplot <- timeplot(
##'   object = fit_cr,
##'   biomarker = response,
##'   id_col = ID,
##'   time_col = time,
##'   time_bin_width = 0.5,
##'   fail_code = 1,
##'   cr_code = 2,
##'   censor_code = 0,
##'   primary_event_label = "Event type 1",
##'   competing_event_label = "Event type 2"
##' )
##'
##' crplot$combined
##' crplot$p_cif
##'
##' }
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##'
##' @export

timeplot <- function(object,
                     biomarker,
                     id_col,
                     time_col,
                     visit_col = NULL,
                     n.obs = 100,
                     seed = 100,
                     window_days = 50,
                     time_bin_width = 180,
                     biomarker_y_lab = NULL,
                     logvar_y_lab = "Log Residual Variance",
                     event_x_lab = NULL,
                     event_y_lab = NULL,
                     x_lab = NULL,
                     x_break_by = NULL,
                     ylim_mean = NULL,
                     ylim_logvar = NULL,
                     event_time_col = NULL,
                     event_status_col = NULL,
                     fail_code = NULL,
                     cr_code = NULL,
                     censor_code = 0,
                     primary_event_label = "Primary event",
                     competing_event_label = NULL,
                     km_title = NULL,
                     cif_title = NULL,
                     center_event_plot = TRUE) {
  
  if (!inherits(object, "jmcs"))
    stop("Use only with 'jmcs' objects.\n")
  
  if (is.null(n.obs) | !is.numeric(n.obs))
    stop("Specify a numbers of subjects to be plotted.")
  
  ydata <- object$ydata
  cdata <- object$cdata
  
  biomarker <- dplyr::enquo(biomarker)
  id_col    <- dplyr::enquo(id_col)
  time_col  <- dplyr::enquo(time_col)
  
  biomarker_nm <- rlang::as_name(biomarker)
  id_nm        <- rlang::as_name(id_col)
  time_nm      <- rlang::as_name(time_col)
  
  if (is.null(x_lab)) {
    x_lab <- time_nm
  }
  
  required_y_cols <- c(biomarker_nm, id_nm, time_nm)
  if (!all(required_y_cols %in% names(ydata))) {
    missing <- setdiff(required_y_cols, names(ydata))
    stop("Missing column(s) in object$ydata: ", paste(missing, collapse = ", "))
  }
  
  visit_nm <- NULL
  
  if (!missing(visit_col) && !is.null(substitute(visit_col))) {
    visit_col <- dplyr::enquo(visit_col)
    visit_nm <- rlang::as_name(visit_col)
    
    if (visit_nm %in% names(ydata)) {
      ydata <- ydata %>%
        dplyr::mutate(.visit_group = .data[[visit_nm]])
      visit_source <- "observed visit"
    } else {
      warning("`visit_col` was supplied but not found in object$ydata: ",
              visit_nm,
              ". Using derived time bins instead.")
      ydata <- ydata %>%
        dplyr::mutate(.visit_group = floor(.data[[time_nm]] / time_bin_width))
      visit_source <- "derived time bin"
      visit_nm <- NULL
    }
  } else {
    ydata <- ydata %>%
      dplyr::mutate(.visit_group = floor(.data[[time_nm]] / time_bin_width))
    visit_source <- "derived time bin"
  }
  
  fitted_vals <- fitted(object, type = "Subject", process = "Longitudinal")
  
  if (length(fitted_vals) != nrow(ydata)) {
    stop("Length of fitted longitudinal values does not match nrow(object$ydata).")
  }
  
  ydata$fitted <- fitted_vals
  
  ydata <- ydata %>%
    dplyr::mutate(
      residual = .data[[biomarker_nm]] - fitted,
      log_resid_sq = log(residual^2 + 1e-8)
    )
  
  traj_res <- plot_biomarker_trajectories(
    data = ydata,
    bio = biomarker_nm,
    id_col = id_nm,
    time_col = time_nm,
    visit_col = if (visit_source == "observed visit") visit_nm else NULL,
    n_ids_mean = n.obs,
    n_ids_bg = n.obs,
    seed = seed,
    time_bin_width = time_bin_width,
    x_lab = x_lab,
    y_lab = biomarker_y_lab,
    x_break_by = x_break_by
  )
  
  p1 <- traj_res$plot
  
  if (!is.null(ylim_mean)) {
    p1 <- p1 + ggplot2::coord_cartesian(ylim = ylim_mean)
  }
  
  subtitle2 <- if (visit_source == "observed visit") {
    "Dark: visit-level variance"
  } else {
    paste0("Dark: time-bin variance (bin width = ", time_bin_width, ")")
  }
  
  ydata_for_p2 <- ydata %>%
    dplyr::filter(.data[[id_nm]] %in% traj_res$ids_mean)
  
  mean_biomarker <- ydata_for_p2 %>%
    dplyr::group_by(.data[[".visit_group"]]) %>%
    dplyr::summarise(
      VisitTime      = mean(.data[[time_nm]], na.rm = TRUE),
      LowerVisitTime = VisitTime - window_days,
      UpperVisitTime = VisitTime + window_days,
      .groups = "drop"
    )
  
  ydataNew <- ydata_for_p2 %>%
    dplyr::left_join(mean_biomarker, by = ".visit_group") %>%
    dplyr::filter(
      .data[[time_nm]] >= LowerVisitTime,
      .data[[time_nm]] <= UpperVisitTime
    )
  
  obs_per_subject_visit <- ydataNew %>%
    dplyr::group_by(.data[[id_nm]], .data[[".visit_group"]]) %>%
    dplyr::summarise(
      n_obs = dplyr::n(),
      .groups = "drop"
    )
  
  mean_biomarkerOneBin <- ydataNew %>%
    dplyr::group_by(.data[[".visit_group"]]) %>%
    dplyr::summarise(
      meanday = mean(.data[[time_nm]], na.rm = TRUE),
      mean_biomarker = mean(.data[[biomarker_nm]], na.rm = TRUE),
      meanres = mean(residual, na.rm = TRUE),
      varres = stats::var(residual, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame()
  
  p2 <- ggplot2::ggplot(
    mean_biomarkerOneBin,
    ggplot2::aes(x = meanday, y = log(varres))
  ) +
    ggplot2::geom_line(linewidth = 1.1, color = "grey5", na.rm = TRUE) +
    ggplot2::geom_point(size = 1.4, color = "grey5", na.rm = TRUE) +
    ggplot2::scale_x_continuous(
      breaks = if (is.null(x_break_by)) {
        ggplot2::waiver()
      } else {
        seq(
          floor(min(mean_biomarkerOneBin$meanday, na.rm = TRUE)),
          ceiling(max(mean_biomarkerOneBin$meanday, na.rm = TRUE)),
          by = x_break_by
        )
      },
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::labs(
      title = "Log residual variance over time",
      subtitle = subtitle2,
      x = x_lab,
      y = logvar_y_lab
    ) +
    theme_timeplot_clean(base_size = 12)
  
  if (!is.null(ylim_logvar)) {
    p2 <- p2 + ggplot2::coord_cartesian(ylim = ylim_logvar)
  }
  
  # Infer event time/status columns from jmcs survival formula when not supplied.
  if (!is.null(object$SurvivalSubmodel)) {
    surv_vars <- all.vars(object$SurvivalSubmodel)
    
    if (is.null(event_time_col)) {
      event_time_col <- surv_vars[1]
    }
    
    if (is.null(event_status_col)) {
      event_status_col <- surv_vars[2]
    }
  }
  
  if (is.null(event_time_col) || is.null(event_status_col)) {
    stop("Could not infer event time/status columns. Please supply `event_time_col` and `event_status_col`.")
  }
  
  if (!all(c(event_time_col, event_status_col) %in% names(cdata))) {
    missing_event_cols <- setdiff(c(event_time_col, event_status_col), names(cdata))
    stop("Missing column(s) in object$cdata: ", paste(missing_event_cols, collapse = ", "))
  }
  
  status_values <- sort(unique(stats::na.omit(cdata[[event_status_col]])))
  event_values <- setdiff(status_values, censor_code)
  
  # Infer fail_code for single-failure models if possible.
  if (!isTRUE(object$CompetingRisk) && is.null(fail_code)) {
    if (length(event_values) == 1) {
      fail_code <- event_values[1]
    } else {
      stop(
        "Could not infer `fail_code`. Non-censor status values are: ",
        paste(event_values, collapse = ", "),
        ". Please supply `fail_code`."
      )
    }
  }
  
  # For competing-risk models, the primary event cannot be inferred from the model alone.
  if (isTRUE(object$CompetingRisk) && is.null(fail_code)) {
    stop(
      "This jmcs object contains competing risks. Please supply `fail_code` to indicate the primary event. ",
      "Observed non-censor status values are: ",
      paste(event_values, collapse = ", ")
    )
  }
  
  # Infer cr_code if there is exactly one other event code.
  if (isTRUE(object$CompetingRisk) && is.null(cr_code)) {
    possible_cr_codes <- setdiff(event_values, fail_code)
    
    if (length(possible_cr_codes) == 1) {
      cr_code <- possible_cr_codes[1]
    } else {
      stop(
        "Could not infer `cr_code`. Non-censor status values are: ",
        paste(event_values, collapse = ", "),
        ". Please supply `cr_code`."
      )
    }
  }
  
  event_res <- NULL
  plot_list <- list(p1, p2)
  
  if (isTRUE(object$CompetingRisk)) {
    event_res <- plot_km_or_cif(
      data = cdata,
      time_col = event_time_col,
      status_col = event_status_col,
      fail_code = fail_code,
      cr_code = cr_code,
      censor_code = censor_code,
      x_lab = if (is.null(event_x_lab)) event_time_col else event_x_lab,
      primary_event_label = primary_event_label,
      competing_event_label = competing_event_label,
      km_title = km_title,
      cif_title = cif_title,
      show_km = FALSE,
      show_cif = TRUE
    )
    
    plot_list <- append(plot_list, list(event_res$cif_plot))
    
  } else {
    event_res <- plot_km_or_cif(
      data = cdata,
      time_col = event_time_col,
      status_col = event_status_col,
      fail_code = fail_code,
      cr_code = NULL,
      censor_code = censor_code,
      x_lab = event_time_col,
      primary_event_label = primary_event_label,
      competing_event_label = NULL,
      km_title = km_title,
      cif_title = NULL,
      show_km = TRUE,
      show_cif = FALSE
    )
    
    plot_list <- append(plot_list, list(event_res$km_plot))
  }
  
  if (!is.null(event_y_lab)) {
    event_res$cif_plot <- event_res$cif_plot + ggplot2::labs(y = event_y_lab)
    event_res$km_plot  <- event_res$km_plot  + ggplot2::labs(y = event_y_lab)
  }
  
  n_plots <- length(plot_list)
  
  if (n_plots == 2) {
    combined <- ggpubr::ggarrange(
      plotlist = plot_list,
      ncol = 2,
      nrow = 1
    )
    
  } else if (n_plots == 3 && isTRUE(center_event_plot)) {
    top_row <- ggpubr::ggarrange(
      plot_list[[1]],
      plot_list[[2]],
      ncol = 2,
      nrow = 1
    )
    
    bottom_row <- ggpubr::ggarrange(
      NULL,
      plot_list[[3]],
      NULL,
      ncol = 3,
      widths = c(0.2, 0.6, 0.2)
    )
    
    combined <- ggpubr::ggarrange(
      top_row,
      bottom_row,
      ncol = 1,
      nrow = 2,
      heights = c(1, 1)
    )
    
  } else if (n_plots == 3) {
    combined <- ggpubr::ggarrange(
      plotlist = plot_list,
      ncol = 1,
      nrow = 3
    )
    
  } else {
    combined <- ggpubr::ggarrange(
      plotlist = plot_list,
      ncol = 2,
      nrow = ceiling(n_plots / 2)
    )
  }
  
  print(combined)
  
  invisible(list(
    combined = combined,
    p1 = p1,
    p2 = p2,
    p_cif = if (!is.null(event_res)) event_res$cif_plot else NULL,
    p_km = if (!is.null(event_res)) event_res$km_plot else NULL,
    cif_fit = if (!is.null(event_res)) event_res$cif_fit else NULL,
    km_fit = if (!is.null(event_res)) event_res$km_fit else NULL,
    event_time_col = event_time_col,
    event_status_col = event_status_col,
    fail_code = fail_code,
    cr_code = cr_code,
    censor_code = censor_code,
    summary_data = mean_biomarkerOneBin,
    obs_per_subject_visit = obs_per_subject_visit,
    grouping_used = visit_source
  ))
}