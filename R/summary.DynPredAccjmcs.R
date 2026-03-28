##' @title Summaries of evaluation metrics for joint models
##' @name summary
##' @aliases summary.DynPredAccjmcs
##' @param object object of class 'DynPredAccjmcs'.
##' @param metric a list to indicate what metric to summarize
##' @param ... Further arguments passed to or from other methods.
##' @return a summary of the list of matrices with conditional probabilities for subjects.
##' @seealso \code{\link{jmcs}, \link{survfitjmcs}}
##' @export
##' 
summary.DynPredAccjmcs <- function(object,
                             metric = c("AUC", "Cindex", "Brier", "MAE", "MAEQ"),
                             digits = 4, ...) {
  if (!inherits(object, "DynPredAccjmcs"))
    stop("Use only with 'DynPredAccjmcs' xs.\n")
  
  metric <- match.arg(metric)
  
  if (!metric %in% object$metrics) {
    stop(paste0("Metric '", metric, "' was not computed in this DynPredAccjmcs object."))
  }
  
  ## helper to check fold success
  .check_cv <- function(x, n.cv) {
    !is.null(x) && length(x) == n.cv && sum(mapply(is.null, x)) == 0
  }
  
  ## helper to average a list of matrices
  .mean_cv_matrix <- function(x, n.cv) {
    total <- 0
    for (j in 1:n.cv) {
      total <- total + x[[j]]
    }
    total / n.cv
  }
  
  if (metric == "AUC") {
    if (!.check_cv(object$AUC.cv, object$n.cv)) {
      stop("The cross validation fails. Please try using a different seed number.")
    }
    
    avg <- .mean_cv_matrix(object$AUC.cv, object$n.cv)
    
    if (object$CompetingRisk) {
      AUC1 <- avg[, 1]
      AUC2 <- avg[, 2]
      ExpectedAUC <- data.frame(object$horizon.time, AUC1, AUC2)
      colnames(ExpectedAUC) <- c("Horizon Time", "AUC1", "AUC2")
    } else {
      AUC <- avg[, 1]
      ExpectedAUC <- data.frame(object$horizon.time, AUC)
      colnames(ExpectedAUC) <- c("Horizon Time", "AUC")
    }
    
    ExpectedAUC[, -1] <- round(ExpectedAUC[, -1], digits)
    cat("\nExpected AUC at the landmark time of", object$landmark.time,
        "\nbased on", object$n.cv, "fold cross validation\n")
    return(ExpectedAUC)
  }
  
  if (metric == "Cindex") {
    if (!.check_cv(object$Cindex.cv, object$n.cv)) {
      stop("The cross validation fails. Please try using a different seed number.")
    }
    
    avg <- .mean_cv_matrix(object$Cindex.cv, object$n.cv)
    
    if (object$CompetingRisk) {
      Cindex1 <- avg[, 1]
      Cindex2 <- avg[, 2]
      ExpectedCindex <- data.frame(object$horizon.time, Cindex1, Cindex2)
      colnames(ExpectedCindex) <- c("Horizon Time", "Cindex1", "Cindex2")
    } else {
      Cindex <- avg[, 1]
      ExpectedCindex <- data.frame(object$horizon.time, Cindex)
      colnames(ExpectedCindex) <- c("Horizon Time", "Cindex")
    }
    
    ExpectedCindex[, -1] <- round(ExpectedCindex[, -1], digits)
    cat("\nExpected Cindex at the landmark time of", object$landmark.time,
        "\nbased on", object$n.cv, "fold cross validation\n")
    return(ExpectedCindex)
  }
  
  if (metric == "Brier") {
    if (!.check_cv(object$Brier.cv, object$n.cv)) {
      stop("The cross validation fails. Please try using a different seed number.")
    }
    
    avg <- .mean_cv_matrix(object$Brier.cv, object$n.cv)
    
    if (object$CompetingRisk) {
      mean.Brier1 <- avg[, 1]
      mean.Brier2 <- avg[, 2]
      ExpectedBrier <- data.frame(object$horizon.time, mean.Brier1, mean.Brier2)
      colnames(ExpectedBrier) <- c("Horizon Time", "Brier Score 1", "Brier Score 2")
    } else {
      mean.Brier <- avg[, 1]
      ExpectedBrier <- data.frame(object$horizon.time, mean.Brier)
      colnames(ExpectedBrier) <- c("Horizon Time", "Brier Score")
    }
    
    ExpectedBrier[, -1] <- round(ExpectedBrier[, -1], digits)
    cat("\nExpected Brier Score at the landmark time of", object$landmark.time,
        "\nbased on", object$n.cv, "fold cross validation\n")
    return(ExpectedBrier)
  }
  
  if (metric == "MAE") {
    if (!.check_cv(object$MAE.cv, object$n.cv)) {
      stop("The cross validation fails. Please try using a different seed number.")
    }
    
    avg <- .mean_cv_matrix(object$MAE.cv, object$n.cv)
    
    if (object$CompetingRisk) {
      mean.MAE1 <- avg[, 1]
      mean.MAE2 <- avg[, 2]
      ExpectedMAE <- data.frame(object$horizon.time, mean.MAE1, mean.MAE2)
      colnames(ExpectedMAE) <- c("Horizon Time", "MAE1", "MAE2")
    } else {
      mean.MAE <- avg[, 1]
      ExpectedMAE <- data.frame(object$horizon.time, mean.MAE)
      colnames(ExpectedMAE) <- c("Horizon Time", "MAE")
    }
    
    ExpectedMAE[, -1] <- round(ExpectedMAE[, -1], digits)
    cat("\nExpected mean absolute error at the landmark time of", object$landmark.time,
        "\nbased on", object$n.cv, "fold cross validation\n")
    return(ExpectedMAE)
  }
  
  if (metric == "MAEQ") {
    if (!.check_cv(object$MAEQ.cv, object$n.cv)) {
      stop("The cross validation fails. Please try using a different seed number.")
    }
    
    if (object$CompetingRisk) {
      out <- as.data.frame(matrix(0, nrow = length(object$horizon.time), ncol = 3))
      out[, 1] <- object$horizon.time
      colnames(out) <- c("Horizon Time", "MAEQ1", "MAEQ2")
      
      for (i in 1:length(object$horizon.time)) {
        for (j in 1:object$n.cv) {
          out[i, 2] <- out[i, 2] + mean(abs(
            object$MAEQ.cv[[j]]$AllCIF1[[i]][, 1] -
              object$MAEQ.cv[[j]]$AllCIF1[[i]][, 2]
          ))
          out[i, 3] <- out[i, 3] + mean(abs(
            object$MAEQ.cv[[j]]$AllCIF2[[i]][, 1] -
              object$MAEQ.cv[[j]]$AllCIF2[[i]][, 2]
          ))
        }
      }
      
      out[, -1] <- out[, -1] / object$n.cv
    } else {
      out <- as.data.frame(matrix(0, nrow = length(object$horizon.time), ncol = 2))
      out[, 1] <- object$horizon.time
      colnames(out) <- c("Horizon Time", "MAEQ")
      
      for (i in 1:length(object$horizon.time)) {
        for (j in 1:object$n.cv) {
          out[i, 2] <- out[i, 2] + mean(abs(
            object$MAEQ.cv[[j]]$AllSurv[[i]][, 1] -
              object$MAEQ.cv[[j]]$AllSurv[[i]][, 2]
          ))
        }
      }
      
      out[, -1] <- out[, -1] / object$n.cv
    }
    
    out[, -1] <- round(out[, -1], digits)
    cat("\nMean absolute error across quantiles of predicted risk scores at the landmark time of",
        object$landmark.time,
        "\nbased on", object$n.cv, "fold cross validation\n")
    return(out)
  }
}