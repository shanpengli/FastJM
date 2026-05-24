##' Simulate joint model data with heterogeneous within-subject variability
##'
##' @title Simulate joint model data with heterogeneous within-subject variability
##' @name simJMWSVdata
##'
##' @description
##' Simulates longitudinal biomarker data with heterogeneous within-subject
##' variability and associated time-to-event data. The survival outcome may be
##' generated either as a single-failure time-to-event outcome or as competing
##' risks data with two failure types.
##'
##' @param seed Integer random seed used for data generation. Default is
##'   \code{100}.
##' @param N Integer specifying the sample size. Default is \code{200}.
##' @param increment Numeric value specifying the increment of visit times for
##'   longitudinal measurements. Default is \code{0.7}.
##' @param beta Numeric vector of true fixed-effect parameters for the
##'   longitudinal mean sub-model.
##' @param tau Numeric vector of true parameters for the within-subject
##'   variability sub-model.
##' @param gamma1 Numeric vector of true fixed-effect parameters in the
##'   cause-specific hazard sub-model for failure type 1.
##' @param gamma2 Numeric vector of true fixed-effect parameters in the
##'   cause-specific hazard sub-model for failure type 2.
##' @param alpha1 Numeric vector of association parameters of individual mean linking the
##'   longitudinal process to the cause-specific hazard for failure type 1.
##' @param alpha2 Numeric vector of association parameters individual mean linking the
##'   longitudinal process to the cause-specific hazard for failure type 2.
##' @param vee1 Numeric association parameter of within-subject
##'   variability linking the longitudinal process to the cause-specific hazard for failure type 1.
##' @param vee2 Numeric association parameter of within-subject
##'   variability linking the longitudinal process to the cause-specific hazard for failure type 2.
##' @param lambda1 Baseline hazard rate for failure type 1. An exponential
##'   baseline hazard with rate \code{lambda1} is assumed.
##' @param lambda2 Baseline hazard rate for failure type 2. An exponential
##'   baseline hazard with rate \code{lambda2} is assumed.
##' @param CL Lower bound of the uniform distribution used to generate censoring
##'   times.
##' @param CU Upper bound of the uniform distribution used to generate censoring
##'   times.
##' @param covbw Variance-covariance matrix for the random effects in the
##'   longitudinal mean and within-subject variability sub-models.
##' @param CR Logical; if \code{TRUE}, competing risks data with two failure
##'   types are simulated. If \code{FALSE}, a single-failure time-to-event
##'   outcome is generated. Default is \code{TRUE}.
##'
##' @return
##' A list with the following components:
##' \describe{
##'   \item{\code{ydatah}}{A long-format data frame containing the simulated
##'   longitudinal biomarker measurements.}
##'   \item{\code{cdatah}}{A data frame containing the simulated event-time data.}
##' }
##'
##' @export

simJMWSVdata <- function(seed = 100, N = 200, increment = 0.7, beta = c(5, 1.5, 2, 1, 2),
                         tau = c(0.5, 0.5, -0.2, 0.2, 0.05),
                         gamma1 = c(1, 0.5, 0.5),
                         gamma2 = c(-0.5, 0.5, 0.25),
                         alpha1 = c(1, 0.7),
                         alpha2 = c(-1, -0.5),
                         vee1 = 0.5,
                         vee2 = -0.5,
                         lambda1 = 0.05,
                         lambda2 = 0.025,
                         CL = 5,
                         CU = 10,
                         covbw = diag(rep(1, 3)),
                         CR = TRUE
                      ) {
  
  set.seed(seed)
  
  bwi <- MASS::mvrnorm(n = N, c(0, 0, 0), covbw, tol = 1e-6, empirical = FALSE)
  
  ##covariate
  Z1 <- sample(c(0, 1), N, replace = TRUE, prob = c(0.5, 0.5))
  Z2 <- runif(N, min = -1, max = 1)
  Z3 <- rnorm(N, mean = 1, sd = 2)
  Z <- cbind(Z1, Z2, Z3)
  
  ## non-informative cencoring time 
  C <- runif(N, min = CL, max = CU)
  #hazard rate of risk1 and risk2
  if (CR) {
    risk1 <- vector()
    risk2 <- vector()
    for (i in 1:N) {
      temp=lambda1*exp(Z[i, ] %*% gamma1 + alpha1%*% bwi[i, 1:2] + vee1*bwi[i, 3])
      risk1[i] <- rexp(1, temp)
      temp=lambda2*exp(Z[i, ] %*% gamma2 + alpha2%*% bwi[i, 1:2] + vee2*bwi[i, 3])
      risk2[i] <- rexp(1, temp)
    }
    survtimeraw <- cbind(risk1, risk2, C)
    
    cmprsk <- vector()
    survtime <- vector()
    for (i in 1:N) {
      if (min(survtimeraw[i, ]) == survtimeraw[i, 1]) {
        cmprsk[i] <- 1
        survtime[i] <- survtimeraw[i, 1]
      } else if (min(survtimeraw[i, ]) == survtimeraw[i, 2]) {
        cmprsk[i] <- 2
        survtime[i] <- survtimeraw[i, 2]
      } else {
        cmprsk[i] <- 0
        survtime[i] <- survtimeraw[i, 3]
      }
    }
    survtimeraw <- as.data.frame(cbind(survtimeraw, survtime, cmprsk))
    table <- as.data.frame(table(survtimeraw$cmprsk)/N*100)
    writeLines(paste0("The censoring rate is: ", table[1, 2], "%"))
    writeLines(paste0("The risk 1 rate is: ", table[2, 2], "%"))
    writeLines(paste0("The risk 2 rate is: ", table[3, 2], "%"))
  } else {
    risk1 <- vector()
    for (i in 1:N) {
      temp=lambda1*exp(Z[i, ] %*% gamma1 + alpha1%*% bwi[i, 1:2] + vee1*bwi[i, 3])
      risk1[i] <- rexp(1, temp)
    }
    survtimeraw <- cbind(risk1, C)
    
    cmprsk <- vector()
    survtime <- vector()
    for (i in 1:N) {
      if (min(survtimeraw[i, ]) == survtimeraw[i, 1]) {
        cmprsk[i] <- 1
        survtime[i] <- survtimeraw[i, 1]
      } else {
        cmprsk[i] <- 0
        survtime[i] <- survtimeraw[i, 2]
      }
    }
    survtimeraw <- as.data.frame(cbind(survtimeraw, survtime, cmprsk))
    table <- as.data.frame(table(survtimeraw$cmprsk)/N*100)
    writeLines(paste0("The censoring rate is: ", table[1, 2], "%"))
    writeLines(paste0("The event rate is: ", table[2, 2], "%"))
  }

  ID <- c(1:N)
  cdata <- cbind(ID, survtimeraw$survtime, survtimeraw$cmprsk, Z)
  colnames(cdata) <- c("ID", "survtime", "cmprsk", "var1", "var2", "var3")
  
  ##fixed effects in longitudinal mean portion
  YdataRaw <- NULL
  for (i in 1:N) {
    ni <- floor(cdata[i, 2]/increment)
    suby <- matrix(0, nrow = ni+1, ncol = 3)
    suby[, 1] <- i
    sd <- sqrt(exp(tau[1] + tau[2]*Z[i, 1] + tau[3]*Z[i, 2] + tau[4]*Z[i, 3] + bwi[i, 3]))
    suby[1, 2] <- beta[1] + beta[2]*Z[i, 1] + beta[3]*Z[i, 2] + beta[4]*Z[i, 3] + bwi[i, 1] + rnorm(1, mean = 0, sd = sd)
    suby[1, 3] <- 0
    if (ni==0) {
      colnames(suby) <- c("ID", "Y", "time")
    } else {
      for (j in 1:ni) {
        sd <- sqrt(exp(tau[1] + tau[2]*Z[i, 1] + tau[3]*Z[i, 2] + tau[4]*Z[i, 3] + tau[5]*j*increment + bwi[i, 3]))
        suby[j+1, 2] <- beta[1] + beta[2]*Z[i, 1] + beta[3]*Z[i, 2] + beta[4]*Z[i, 3] + beta[5]*j*increment + 
          bwi[i, 1] + bwi[i, 2]*j + rnorm(1, mean = 0, sd = sd) 
        suby[j+1, 3] <- j
      }
    }
    YdataRaw <- rbind(YdataRaw, suby)
    
  }
  Z <- cbind(ID, Z)
  Z <- as.data.frame(Z)
  YdataRaw <- as.data.frame(YdataRaw)
  ydata <- dplyr::left_join(YdataRaw, Z, by = "ID")
  cdata <- as.data.frame(cdata)
  a <- list(cdata, ydata)
  names(a) <- c("cdatah", "ydatah")
  return(a)
}




