plot_biomarker_trajectories <- function(data,
                                        bio,
                                        id_col = "ID",
                                        time_col,
                                        visit_col = NULL,
                                        time_bin_width = 180,
                                        n_ids_mean = 100,
                                        plot_title = NULL,
                                        y_lab = NULL,
                                        x_lab = time_col,
                                        n_ids_bg = 30,
                                        seed = 123,
                                        x_break_by = NULL,
                                        show_ribbon = TRUE,
                                        ribbon_alpha = 0.15,
                                        line_alpha = 0.30,
                                        line_width = 0.3,
                                        mean_width = 1.1,
                                        mean_point_size = 1.4,
                                        bg_color = "grey35",
                                        mean_color = "grey5") {
  stopifnot(is.data.frame(data))
  
  required_cols <- c(id_col, time_col, bio)
  if (!all(required_cols %in% names(data))) {
    missing <- setdiff(required_cols, names(data))
    stop("Missing column(s) in `data`: ", paste(missing, collapse = ", "))
  }
  
  if (!is.null(visit_col) && !visit_col %in% names(data)) {
    warning("`visit_col` was supplied but not found in `data`: ",
            visit_col,
            ". Using derived time bins instead.")
    visit_col <- NULL
  }
  
  ID   <- id_col
  TIME <- time_col
  BIO  <- bio
  
  if (!is.null(seed)) set.seed(seed)
  
  temp0 <- data
  
  if (is.null(visit_col)) {
    temp0 <- temp0 %>%
      dplyr::mutate(.visit_group = floor(.data[[TIME]] / time_bin_width))
    VISIT <- ".visit_group"
  } else {
    VISIT <- visit_col
  }
  
  id_tbl <- temp0 %>% dplyr::distinct(.data[[ID]])
  n_take_mean <- min(n_ids_mean, nrow(id_tbl))
  
  ids_mean <- id_tbl %>%
    dplyr::slice_sample(n = n_take_mean) %>%
    dplyr::pull(.data[[ID]])
  
  temp <- temp0 %>%
    dplyr::filter(.data[[ID]] %in% ids_mean) %>%
    dplyr::arrange(.data[[ID]], .data[[TIME]]) %>%
    dplyr::group_by(.data[[ID]]) %>%
    #dplyr::filter(sum(!is.na(.data[[BIO]])) >= 2) %>%
    dplyr::ungroup()
  
  id_tbl_bg <- temp %>% dplyr::distinct(.data[[ID]])
  n_take_bg <- min(n_ids_bg, nrow(id_tbl_bg))
  
  ids_bg <- id_tbl_bg %>%
    dplyr::slice_sample(n = n_take_bg) %>%
    dplyr::pull(.data[[ID]])
  
  temp_bg <- temp %>%
    dplyr::filter(.data[[ID]] %in% ids_bg)
  
  sum_syn <- temp %>%
    dplyr::group_by(.data[[VISIT]]) %>%
    dplyr::summarise(
      mean_day   = mean(.data[[TIME]], na.rm = TRUE),
      mean_value = mean(.data[[BIO]], na.rm = TRUE),
      se         = stats::sd(.data[[BIO]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[BIO]]))),
      lo         = mean_value - 1.96 * se,
      hi         = mean_value + 1.96 * se,
      n          = sum(!is.na(.data[[BIO]])),
      .groups    = "drop"
    )
  
  x_max <- max(temp[[TIME]], na.rm = TRUE)
  x_min <- min(temp[[TIME]], na.rm = TRUE)
  breaks_x <- if (is.null(x_break_by)) {
    ggplot2::waiver()
  } else {
    seq(floor(x_min), ceiling(x_max), by = x_break_by)
  }
  
  subtitle_txt <- if (!is.null(visit_col)) {
    "Dark line: visit-level mean, grey lines: individual patients"
  } else {
    paste0("Dark: time-bin mean (bin width = ", time_bin_width,
           "), grey lines: individual patients")
  }
  
  p <- ggplot2::ggplot()
  
  if (isTRUE(show_ribbon)) {
    p <- p + ggplot2::geom_ribbon(
      data = sum_syn,
      ggplot2::aes(x = mean_day, ymin = lo, ymax = hi),
      alpha = ribbon_alpha
    )
  }
  
  p <- p +
    ggplot2::geom_path(
      data = temp_bg,
      ggplot2::aes(x = .data[[TIME]], y = .data[[BIO]], group = .data[[ID]]),
      linewidth = line_width,
      alpha = line_alpha,
      color = bg_color,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      data = sum_syn,
      ggplot2::aes(x = mean_day, y = mean_value),
      linewidth = mean_width,
      color = mean_color
    ) +
    ggplot2::geom_point(
      data = sum_syn,
      ggplot2::aes(x = mean_day, y = mean_value),
      size = mean_point_size,
      color = mean_color
    ) +
    ggplot2::scale_x_continuous(
      breaks = breaks_x,
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::labs(
      title = if (is.null(plot_title)) {
        # paste0("Trajectories of ", BIO, " (n=", length(unique(temp[[ID]])),
        #        " IDs; ", length(unique(temp_bg[[ID]])), " trajectories)")
        paste0("Trajectories of ", BIO, " (n=", length(unique(temp[[ID]])), ")")
      } else {
        plot_title
      },
      subtitle = subtitle_txt,
      x = x_lab,
      y = BIO
    ) +
    theme_timeplot_clean(base_size = 12)
  
  list(
    plot = p,
    ids_mean = ids_mean,
    ids_bg = ids_bg,
    summary_data = sum_syn,
    data_mean = temp,
    data_bg = temp_bg
  )
}


########## KM and CIF helper function

plot_km_or_cif <- function(data,
                           time_col,
                           status_col,
                           fail_code,
                           cr_code = NULL,
                           censor_code = 0,
                           x_lab = NULL,
                           primary_event_label = "Primary event",
                           competing_event_label = NULL,
                           km_title = NULL,
                           cif_title = NULL,
                           show_km = FALSE,
                           show_cif = FALSE) {
  stopifnot(is.data.frame(data))
  
  if (is.null(x_lab)) {
    x_lab <- time_col
  }
  
  if (is.null(fail_code)) {
    stop("Please supply `fail_code`.")
  }
  
  if (isTRUE(show_cif) && is.null(cr_code)) {
    stop("`cr_code` is required when `show_cif = TRUE`.")
  }
  
  if (isTRUE(show_cif) && is.null(competing_event_label)) {
    competing_event_label <- "Competing event"
  }
  
  required_cols <- c(time_col, status_col)
  if (!all(required_cols %in% names(data))) {
    missing <- setdiff(required_cols, names(data))
    stop("Missing column(s) in event data: ", paste(missing, collapse = ", "))
  }
  
  d <- data[, c(time_col, status_col), drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  names(d) <- c("time", "status")
  
  allowed_status <- c(censor_code, fail_code)
  if (!is.null(cr_code)) {
    allowed_status <- c(allowed_status, cr_code)
  }
  
  bad_status <- setdiff(unique(d$status), allowed_status)
  if (length(bad_status) > 0) {
    stop(
      "`status` contains values outside the supplied event/censoring codes: ",
      paste(bad_status, collapse = ", "),
      ". Supplied codes were: ",
      paste(allowed_status, collapse = ", ")
    )
  }
  
  out <- list(
    cif_fit = NULL,
    cif_plot = NULL,
    km_fit = NULL,
    km_plot = NULL
  )
  
  if (isTRUE(show_cif)) {
    d$cif_status <- factor(
      d$status,
      levels = c(censor_code, fail_code, cr_code),
      labels = c("Censored", primary_event_label, competing_event_label)
    )
    
    cif_fit <- tidycmprsk::cuminc(
      survival::Surv(time, cif_status) ~ 1,
      data = d
    )
    
    p_cif <-
      ggsurvfit::ggcuminc(
        cif_fit,
        outcome = c(primary_event_label, competing_event_label)
      ) +
      ggsurvfit::scale_ggsurvfit() +
      ggplot2::labs(
        title = if (is.null(cif_title)) {
          paste("Cumulative incidence of", primary_event_label, "and", competing_event_label)
        } else {
          cif_title
        },
        x = x_lab,
        y = "Cumulative incidence rate"
      ) +
      theme_timeplot_clean(base_size = 12)
    
    out$cif_fit <- cif_fit
    out$cif_plot <- p_cif
  }
  
  if (isTRUE(show_km)) {
    d$km_status <- as.integer(d$status == fail_code)
    
    km_fit <- survival::survfit(
      survival::Surv(time, km_status) ~ 1,
      data = d
    )
    
    p_km <-
      ggsurvfit::ggsurvfit(km_fit) +
      ggsurvfit::add_confidence_interval() +
      ggsurvfit::add_censor_mark() +
      ggsurvfit::scale_ggsurvfit() +
      ggplot2::labs(
        title = if (is.null(km_title)) {
          paste0("Kaplan-Meier for ", primary_event_label)
        } else {
          km_title
        },
        x = x_lab,
        y = "Survival probability"
      ) +
      theme_timeplot_clean(base_size = 12)
    
    out$km_fit <- km_fit
    out$km_plot <- p_km
  }
  
  out
}

theme_timeplot_clean <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.7),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 4)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 4)),
      
      plot.title = ggplot2::element_text(face = "bold", color = "black"),
      plot.subtitle = ggplot2::element_text(color = "grey30")
    )
}

format_runtime <- function(x, digits = 2) {
  sec <- as.numeric(x, units = "secs")
  if (is.na(sec)) return(NA_character_)
  if (sec < 60) {
    paste0(round(sec, digits), " seconds")
  } else if (sec < 3600) {
    paste0(round(sec / 60, digits), " minutes")
  } else {
    paste0(round(sec / 3600, digits), " hours")
  }
}
