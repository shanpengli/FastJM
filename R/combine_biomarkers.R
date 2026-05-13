##' @title Combine biomarker measurements across multiple data frames
##' @name combine_biomarkers
##' @aliases combine_biomarkers
##' @description
##' Combines specified biomarker measurements across a list of data frames and
##' returns a unified long-format data set restricted to subjects with at least
##' one recorded measurement for EVERY requested biomarker.
##'
##' The function searches each provided data frame for the biomarker names listed
##' in \code{biomarkers}, extracts the matching measurement rows, and standardizes
##' the output into a common structure. Subjects are retained only if they have
##' at least one non-missing measurement for every requested biomarker across the
##' supplied data frames. It is intended for settings where biomarker measurements 
##' may be distributed across multiple data frames, and where each data frame may 
##' contain one or more of the requested biomarkers. The returned data are in long 
##' format so that repeated measurements per subject and per biomarker are preserved. 
##' No collapsing, averaging, or other within-subject summarization is performed.
##'
##' @param data_list A non-empty list of data frames containing subject-level
##'   biomarker measurements.
##' @param biomarkers A character vector of biomarker column names to search for
##'   across the supplied data frames.
##' @param id_col A character string giving the subject ID column name. This
##'   column must be present in every data frame in \code{data_list}.
##' @param time_col An optional character string giving the time variable column
##'   name, such as days from baseline. If supplied, this column must be present
##'   in every data frame in \code{data_list}, and it is standardized to
##'   \code{"time"} in the returned long-format output.
##' @param dataset_names An optional character vector of names for the supplied
##'   data frames. If omitted, the function uses \code{names(data_list)} when
##'   available; otherwise it creates default names of the form
##'   \code{"dataset_1"}, \code{"dataset_2"}, and so on.
##'
##' @return A list with the following components:
##' \describe{
##'   \item{\code{combined_long}}{A long-format data frame containing only
##'   included subjects. The output includes the subject ID column, a
##'   \code{dataset} column identifying the source data frame, a
##'   \code{biomarker} column identifying the biomarker name, a standardized
##'   \code{value} column containing the measurement values, and a standardized
##'   \code{time} column if \code{time_col} was supplied.}
##'   \item{\code{included_ids}}{A vector of subject IDs retained in the final
##'   combined data. These are the subjects with at least one non-missing
##'   measurement for every requested biomarker.}
##'   \item{\code{ids_by_biomarker}}{A named list giving the subject IDs with at
##'   least one non-missing measurement for each biomarker.}
##'   \item{\code{presence_summary}}{A data frame summarizing biomarker presence
##'   among the included subjects, with one row per subject and one logical
##'   column per biomarker.}
##' }
##'
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##'
##' @examples
##' df_a <- data.frame(
##'   ID = c(1, 1, 2, 3, 4, 5),
##'   day = c(0, 30, 0, 0, 0, 0),
##'   sbp = c(120, 125, 130, NA, 110, 118)
##' )
##'
##' df_b <- data.frame(
##'   ID = c(1, 2, 2, 3, 5, 6),
##'   day = c(0, 0, 20, 0, 0, 0),
##'   dbp = c(80, 85, 84, 90, NA, 88)
##' )
##'
##' df_c <- data.frame(
##'   ID = c(1, 1, 2, 4, 5, 5),
##'   day = c(0, 40, 0, 0, 0, 10),
##'   bpm = c(70, 72, 68, 75, 77, 79)
##' )
##'
##' res <- combine_biomarkers(
##'   data_list = list(df_a, df_b, df_c),
##'   biomarkers = c("sbp", "dbp", "bpm"),
##'   id_col = "ID",
##'   time_col = "day"
##' )
##'
##' res$included_ids
##' head(res$combined_long)
##' res$presence_summary
##'
##' @export
combine_biomarkers <- function(data_list,
                                     biomarkers,
                                     id_col,
                                     time_col = NULL,
                                     dataset_names = NULL) {
  # Stop conditions
  if (!is.list(data_list) || length(data_list) == 0) {
    stop("`data_list` must be a non-empty list of data frames.")
  }
  if (length(biomarkers) == 0) {
    stop("`biomarkers` must contain at least one biomarker name.")
  }
  if (is.null(dataset_names)) {
    dataset_names <- names(data_list)
    if (is.null(dataset_names) || any(dataset_names == "")) {
      dataset_names <- paste0("dataset_", seq_along(data_list))
    }
  }
  if (length(dataset_names) != length(data_list)) {
    stop("`dataset_names` must have same length as `data_list`.")
  }
  
  # Checking ID/time columns
  for (i in seq_along(data_list)) {
    dat <- data_list[[i]]
    if (!is.data.frame(dat)) {
      stop("Element ", i, " of `data_list` is not a data frame.")
    }
    if (!id_col %in% names(dat)) {
      stop("`id_col` not found in dataset: ", dataset_names[i])
    }
    if (!is.null(time_col) && !time_col %in% names(dat)) {
      stop("`time_col` not found in dataset: ", dataset_names[i])
    }
  }
  
  # Gathering all requested biomarker rows from ALL datasets
  long_list <- list()
  for (i in seq_along(data_list)) {
    dat <- data_list[[i]]
    dat_name <- dataset_names[i]
    
    bios_here <- intersect(biomarkers, names(dat))
    if (length(bios_here) == 0) next
    
    for (bio in bios_here) {
      keep_cols <- c(id_col, bio)
      if (!is.null(time_col)) keep_cols <- c(keep_cols, time_col)
      
      tmp <- dat[, keep_cols, drop = FALSE]
      names(tmp)[names(tmp) == bio] <- "value"
      if (!is.null(time_col)) {
        names(tmp)[names(tmp) == time_col] <- "time"
      }
      
      tmp$dataset <- dat_name
      tmp$biomarker <- bio
      
      if (is.null(time_col)) {
        tmp <- tmp[, c(id_col, "dataset", "biomarker", "value"), drop = FALSE]
      } else {
        tmp <- tmp[, c(id_col, "dataset", "biomarker", "time", "value"), drop = FALSE]
      }
      
      long_list[[length(long_list) + 1]] <- tmp
    }
  }
  
  if (length(long_list) == 0) {
    stop("None of the requested biomarkers were found in the supplied data frames.")
  }
  
  combined_long <- do.call(rbind, long_list)
  rownames(combined_long) <- NULL
  
  # Subjects with >=1 non-missing value for each biomarker
  ids_by_biomarker <- lapply(biomarkers, function(bio) {
    unique(combined_long[[id_col]][
      combined_long$biomarker == bio & !is.na(combined_long$value)
    ])
  })
  names(ids_by_biomarker) <- biomarkers
  
  included_ids <- Reduce(intersect, ids_by_biomarker)
  
  combined_long_filtered <- combined_long[
    combined_long[[id_col]] %in% included_ids,
    ,
    drop = FALSE
  ]
  
  # Presence summary
  presence_summary <- data.frame(
    ID_tmp = included_ids,
    stringsAsFactors = FALSE
  )
  names(presence_summary)[1] <- id_col
  
  for (bio in biomarkers) {
    presence_summary[[bio]] <- presence_summary[[id_col]] %in% ids_by_biomarker[[bio]]
  }
  
  list(
    combined_long = combined_long_filtered,
    included_ids = included_ids,
    ids_by_biomarker = ids_by_biomarker,
    presence_summary = presence_summary
  )
}