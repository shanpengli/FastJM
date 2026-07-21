##' Simulate multivariate longitudinal and time-to-event data from a joint model
##'
##' @title Data simulation for multivariate joint models conditional on a landmark time
##' @name simmvJMdatalm
##'
##' @description
##' Simulate data from a joint model with multiple longitudinal biomarkers and
##' either a single time-to-event outcome or competing risks outcomes. The
##' longitudinal biomarkers are generated from biomarker-specific linear
##' mixed-effects models, allowing different random-effects structures across
##' biomarkers. The event-time process is generated from cause-specific
##' proportional hazards models with associations induced through the
##' biomarker-specific random-effects contribution or current biomarker value at
##' a landmark time.
##'
##' @param seed Integer specifying the random seed used for data simulation.
##' @param N Integer specifying the sample size.
##' @param increment Numeric scalar specifying the time increment between
##' scheduled longitudinal measurements.
##' @param beta A list of numeric vectors specifying the true fixed-effect
##' parameters for the longitudinal submodels. Each list component corresponds
##' to one biomarker. The entries are interpreted as the intercept, covariate
##' effects, and time effect; if quadratic time is included, the final entry
##' corresponds to the quadratic time effect.
##' @param sigma Numeric vector specifying the measurement error variances for
##' the longitudinal biomarkers. The length of \code{sigma} should equal the
##' number of biomarkers.
##' @param gamma1 Numeric vector specifying the true survival fixed-effect
##' parameters for failure type 1.
##' @param gamma2 Numeric vector specifying the true survival fixed-effect
##' parameters for failure type 2. Used when \code{CR = TRUE}.
##' @param alpha1 A list of numeric vectors specifying the true association
##' parameters between the longitudinal biomarkers and failure type 1. Each
##' component corresponds to one biomarker and should have length matching the
##' corresponding entry of \code{pRE}.
##' @param alpha2 A list of numeric vectors specifying the true association
##' parameters between the longitudinal biomarkers and failure type 2. Each
##' component corresponds to one biomarker and should have length matching the
##' corresponding entry of \code{pRE}. Used when \code{CR = TRUE}.
##' @param lambda1 Numeric scalar specifying the constant baseline hazard rate
##' for failure type 1. An exponential baseline hazard is assumed.
##' @param lambda2 Numeric scalar specifying the constant baseline hazard rate
##' for failure type 2. An exponential baseline hazard is assumed. Used when
##' \code{CR = TRUE}.
##' @param CL Numeric scalar specifying the lower bound of the uniform
##' distribution used to generate censoring times.
##' @param CU Numeric scalar specifying the upper bound of the uniform
##' distribution used to generate censoring times.
##' @param covb Variance-covariance matrix for the subject-specific random
##' effects. Its dimension must be equal to \code{sum(pRE)}.
##' @param missprob Numeric scalar between 0 and 1 specifying the probability
##' that a scheduled longitudinal observation is missing. Default is \code{0}.
##' @param landmark Logical; if \code{TRUE}, event times and censoring times are
##' generated conditional on survival beyond a landmark time \code{s}. Default
##' is \code{TRUE}.
##' @param s Numeric scalar specifying the landmark time. Required when
##' \code{landmark = TRUE}. Default is \code{NULL}.
##' @param method Character string specifying how the longitudinal process is
##' linked to the event-time process. Options are \code{"presentlp"} and
##' \code{"present"}. If \code{method = "presentlp"}, the event-time model
##' depends on the subject-specific random-effects contribution at the landmark
##' time. If \code{method = "present"}, the event-time model depends on the full
##' current biomarker value at the landmark time, including both fixed- and
##' random-effects components. Default is \code{"presentlp"}.
##' @param pRE Integer vector specifying the number of random effects for each
##' biomarker. For example, \code{pRE = c(1, 2)} specifies a random-intercept
##' model for biomarker 1, a random intercept-and-slope model for biomarker 2.
##' If \code{NULL}, a random-intercept model is assumed for
##' each biomarker. The sum of \code{pRE} must equal the dimension of
##' \code{covb}. Default is \code{NULL}.
##' @param CR Logical; if \code{TRUE}, simulate competing risks data with two
##' failure types. If \code{FALSE}, simulate a single failure type with
##' independent censoring. Default is \code{TRUE}.
##'
##' @return
##' A list with two elements:
##' \item{mvcdata}{A data frame containing survival data, including subject ID,
##' observed event or censoring time, event indicator, and baseline covariates.}
##' \item{mvydata}{A long-format data frame containing longitudinal biomarker
##' measurements, visit times, subject ID, and baseline covariates.}
##'
##' @details
##' The function generates two baseline covariates, \code{X1} and \code{X2}.
##' The longitudinal biomarkers are generated from biomarker-specific linear
##' mixed-effects models. The number of random effects for each biomarker is
##' controlled by \code{pRE}. Currently, \code{pRE} values of 1 and 2 are
##' supported, corresponding respectively to random intercept, random intercept
##' and slope.
##'
##' When \code{landmark = TRUE}, event and censoring times are generated after
##' the landmark time \code{s}. Under \code{method = "presentlp"}, the survival
##' model uses only the random-effects contribution evaluated at \code{s}. Under
##' \code{method = "present"}, the survival model uses the full current
##' biomarker value at \code{s}.
##'
##' @examples
##' \dontrun{
##' dat <- simmvJMdatalm(
##'   seed = 100,
##'   N = 5000,
##'   increment = 0.7,
##'   beta = list(
##'     beta1 = c(5, 1.5, 2, 1),
##'     beta2 = c(10, 1, 2, 1),
##'     beta3 = c(8, 1.2, 1.5, 0.8)
##'   ),
##'   sigma = rep(1, 3),
##'   gamma1 = c(1, 0.5),
##'   gamma2 = c(-0.5, 0.5),
##'   alpha1 = list(
##'     alpha11 = -0.5,
##'     alpha12 = c(0.5, 0.7),
##'     alpha13 = c(0.3, 0.4)
##'   ),
##'   alpha2 = list(
##'     alpha21 = 0.5,
##'     alpha22 = c(0.5, 0.8),
##'     alpha23 = c(0.3, 0.1)
##'   ),
##'   lambda1 = 0.05,
##'   lambda2 = 0.05,
##'   covb = diag(c(5, 10, 1, 10, 1)),
##'   pRE = c(1, 2, 2),
##'   s = 2,
##'   landmark = TRUE,
##'   CR = TRUE
##' )
##'
##' mvydata <- dat$mvydata
##' mvcdata <- dat$mvcdata
##' }
##'
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @export

simmvJMdatalm <- function(seed = 100, N = 200, increment = 0.7, beta = list(beta1 = c(5, 1.5, 2, 1),
                                                                          beta2 = c(10, 1, 2, 1)),
                        sigma = c(0.5, 0.5),
                        gamma1 = c(1, 0.5),
                        gamma2 = c(-0.5, 0.5),
                        alpha1 = list(alpha11 = c(0.5),
                                      alpha12 = c(-0.5)),
                        alpha2 = list(alpha21 = c(0.5),
                                      alpha22 = c(-0.5)),
                        lambda1 = 0.05,
                        lambda2 = 0.025,
                        CL = 5,
                        CU = 10,
                        covb = diag(rep(1, 4)),
                        missprob = 0, landmark = TRUE,
                        s = NULL, method = "presentlp", pRE = NULL, CR = TRUE){
  

  
  set.seed(seed)
  
  # for(g in 1:length(beta)){
  #   # assume all random effects have same length for now; will change later
  #   pRE[g] <- length(covb)/length(beta)
  #   pREtotal = pREtotal + pRE[g]
  # }
  
  if (is.null(pRE)) {
    pRE <- rep(1, length(beta))  # assume random intercept only
  }
  
  if (length(pRE) != length(beta)) {
    stop("pRE must have one entry per biomarker.")
  }
  
  if (sum(pRE) != ncol(covb)) {
    stop("sum(pRE) must equal ncol(covb).")
  }
  index = 0
  
  # for(g in 1:length(beta)){
  #   # assume all random effects have same length for now; will change later
  #   # pRE[g] <- length(covb)/length(beta)
  #   pREtotal = pREtotal + pRE[g]
  # }
  pREtotal <- sum(pRE)
  
  
  if(CR == TRUE){
    
    bi <- MASS::mvrnorm(n = N, rep(0, pREtotal), covb, tol = 1e-6, empirical = FALSE)
    
    ##covariate
    X1 <- sample(c(0, 1), N, replace = TRUE, prob = c(0.5, 0.5))
    X2 <- runif(N, min = -5, max = 5)
    X <- cbind(X1, X2)
    
    #hazard rate of risk1 and risk2
    ## non-informative censoring time 
    if (landmark) {
      if (is.null(s)) stop("Must specify landmark time s")
      C <- s + runif(N, min = CL, max = CU)
    } else {
      C <- runif(N, min = CL, max = CU)
    }
    risk1 <- vector()
    risk2 <- vector()
    for (i in 1:N) {
      
      if (landmark && is.null(s)) stop("Must specify landmark time s")
      
      s_eval <- if (landmark) s else 0
      sum.alpha1i <- 0
      sum.alpha2i <- 0
      index <- 0
      
      for (g in 1:length(beta)) {
        
        subbeta <- beta[[g]]
        alpha1g <- alpha1[[g]]
        alpha2g <- alpha2[[g]]
        bi_g <- matrix(bi[i, (index + 1):(index + pRE[g])], ncol = 1)
        
        Xi_s <- matrix(c(1, X[i, ], s_eval), nrow = 1)
        
        if (length(subbeta) != ncol(Xi_s)) {
          stop("Length of beta must match intercept + covariates + time.")
        }
        
        if (pRE[g] == 1) {
          Zi_s <- matrix(1, nrow = 1)
        } else if (pRE[g] == 2) {
          Zi_s <- matrix(c(1, s_eval), nrow = 1)
        } else if (pRE[g] == 3) {
          Zi_s <- matrix(c(1, s_eval, s_eval^2), nrow = 1)
        } else {
          stop("Unsupported pRE[g]")
        }
        
        if (method == "presentlp") {
          m_is <- as.numeric(Zi_s %*% bi_g)
          
        } else if (method == "present") {
          m_is <- as.numeric(Xi_s %*% subbeta + Zi_s %*% bi_g)
          
        } else {
          stop("method must be 'presentlp' or 'present'")
        }
        sum.alpha1i <- sum.alpha1i + alpha1g * m_is
        sum.alpha2i <- sum.alpha2i + alpha2g * m_is
        index <- index + pRE[g]
      }
      
      temp <- lambda1 * exp(X[i, ] %*% gamma1 + sum.alpha1i)
      risk1[i] <- if (landmark) s + rexp(1, rate = temp) else rexp(1, rate = temp)
      temp <- lambda2 * exp(X[i, ] %*% gamma2 + sum.alpha2i)
      risk2[i] <- if (landmark) s + rexp(1, rate = temp) else rexp(1, rate = temp)
      # add landmark time s
      
    }
  
  
    # adjust censoring consider uniform for c
    # look at rsik 1 and risk2 30-50% approx
    
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
    ID <- c(1:N)
    cdata <- cbind(ID, survtimeraw$survtime, survtimeraw$cmprsk, X)

    
    colnames(cdata) <- c("ID", "survtime", "cmprsk", "X21", "X22")
    
    ##fixed effects in longitudinal mean portion
    Ydata <- data.frame(c(1:N))
    colnames(Ydata) <- "ID"
    
    
    index = 0
    
    for (g in 1:length(beta)) {
      
      YdataRaw <- NULL
      subbeta <- beta[[g]]
      Z <- 1
      for (i in 1:N) {
        
        ni <- floor(cdata[i, 2]/increment)
        Visit <- sample(c(0, 1), ni, replace = TRUE, prob = c(missprob, 1 - missprob))
        suby <- matrix(0, nrow = ni+1, ncol = 3)
        suby[, 1] <- i
        
        Xi_0 <- matrix(c(1, X[i, ], 0), nrow = 1)
        
        if (pRE[g] == 1) {
          Z0 <- matrix(1, nrow = 1)
        } else if (pRE[g] == 2) {
          Z0 <- matrix(c(1, 0), nrow = 1)
        } else if (pRE[g] == 3) {
          Z0 <- matrix(c(1, 0, 0^2), nrow = 1)
        } else {
          stop("Unsupported pRE[g]")
        }
        
        suby[1, 3] <- as.numeric(Xi_0 %*% subbeta) +
          as.numeric(Z0 %*% bi[i, (index + 1):(index + pRE[g])]) +
          rnorm(1, mean = 0, sd = sqrt(sigma[g]))
        
        suby[1, 2] <- 0
        
        if (ni==0) {
          
          colnames(suby) <- c("ID", "time", "Y")
          colnames(suby)[3] <- paste0("Y", g) 
          
        } else {
          
          for (j in 1:ni) {
            
            if (pRE[g] == 1) {
              Z <- 1
            } else if (pRE[g] == 2) {
              Z <- c(1, j*increment)
            } else if (pRE[g] == 3) {
              Z <- c(1, j*increment, (j*increment)^2)
            }
            
            if (Visit[j] == 1) {
              
              time_ij <- j * increment
              Xi_t <- matrix(c(1, X[i, ], time_ij), nrow = 1)
              
              suby[j+1, 3] <- as.numeric(Xi_t %*% subbeta) +
                as.numeric(Z %*% bi[i, (index + 1):(index + pRE[g])]) +
                rnorm(1, mean = 0, sd = sqrt(sigma[g]))
              
            } else {
              suby[j+1, 3] <- NA
            }
            
            suby[j+1, 2] <- j*increment
            
          }
          
          colnames(suby) <- c("ID", "time", "Y")
          colnames(suby)[3] <- paste0("Y", g) 
          
        }
        YdataRaw <- rbind(YdataRaw, suby)
      }
      YdataRaw <- as.data.frame(YdataRaw)
      
      if(g == 1) {
        
        Ydata <- dplyr::full_join(Ydata, YdataRaw, by = "ID")
      } else {
        Ydata <- dplyr::full_join(Ydata, YdataRaw, by = c("ID", "time"))
      }
      
      index = index + pRE[g]
      
    }
    
    X <- cbind(ID, X)
    X <- as.data.frame(X)
    colnames(X)[2:3] <- c("X11", "X12")
    Ydata <- as.data.frame(Ydata)
    ydata <- dplyr::left_join(Ydata, X, by = "ID")
    cdata <- as.data.frame(cdata)
    
    a <- list(cdata, ydata)
    names(a) <- c("mvcdata", "mvydata")
    return(a)
    
  } else{
    
    

    bi <- MASS::mvrnorm(n = N, rep(0, pREtotal), covb, tol = 1e-6, empirical = FALSE)
   
    
    ##covariate
    X1 <- sample(c(0, 1), N, replace = TRUE, prob = c(0.5, 0.5))
    X2 <- runif(N, min = -5, max = 5)
    X <- cbind(X1, X2)
    
    #hazard rate of risk1 and risk2
    ## non-informative censoring time 
    if (landmark) {
      if (is.null(s)) stop("Must specify landmark time s")
      C <- s + runif(N, min = CL, max = CU)
    } else {
      C <- runif(N, min = CL, max = CU)
    }
    risk1 <- vector()
    
    for (i in 1:N) {
      
      if (landmark && is.null(s)) stop("Must specify landmark time s")
      
      s_eval <- if (landmark) s else 0
      sum.alpha1i <- 0
      index <- 0
      
      for (g in 1:length(beta)) {
        
        subbeta <- beta[[g]]
        alpha1g <- alpha1[[g]]
        bi_g <- matrix(bi[i, (index + 1):(index + pRE[g])], ncol = 1)
        
        Xi_s <- matrix(c(1, X[i, ], s_eval), nrow = 1)
        
        if (length(subbeta) != ncol(Xi_s)) {
          stop("Length of beta must match intercept + covariates + time.")
        }
        
        if (pRE[g] == 1) {
          Zi_s <- matrix(1, nrow = 1)
        } else if (pRE[g] == 2) {
          Zi_s <- matrix(c(1, s_eval), nrow = 1)
        } else if (pRE[g] == 3) {
          Zi_s <- matrix(c(1, s_eval, s_eval^2), nrow = 1)
        } else {
          stop("Unsupported pRE value")
        }
        
        if (method == "presentlp") {
          m_i <- as.numeric(Zi_s %*% bi_g)
          
        } else if (method == "present") {
          m_i <- as.numeric(Xi_s %*% subbeta + Zi_s %*% bi_g)
          
        } else {
          stop("method must be 'presentlp' or 'present'")
        }
        sum.alpha1i <- sum.alpha1i + alpha1g * m_i
        index <- index + pRE[g]
      }
      
      temp <- lambda1 * exp(X[i, ] %*% gamma1 + sum.alpha1i)
      risk1[i] <- if (landmark) s + rexp(1, rate = temp) else rexp(1, rate = temp)
  
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
    writeLines(paste0("The risk 1 rate is: ", table[2, 2], "%"))
    ID <- c(1:N)
    cdata <- cbind(ID, survtimeraw$survtime, survtimeraw$cmprsk, X)

    colnames(cdata) <- c("ID", "survtime", "cmprsk", "X21", "X22")
    
    ##fixed effects in longitudinal mean portion
    Ydata <- data.frame(c(1:N))
    colnames(Ydata) <- "ID"
    index = 0
    for (g in 1:length(beta)) {
      
      YdataRaw <- NULL
      subbeta <- beta[[g]]
      Z <- 1
      for (i in 1:N) {
        
        ni <- floor(cdata[i, 2]/increment)
        Visit <- sample(c(0, 1), ni, replace = TRUE, prob = c(missprob, 1 - missprob))
        suby <- matrix(0, nrow = ni+1, ncol = 3)
        suby[, 1] <- i
        
        Xi_0 <- matrix(c(1, X[i, ], 0), nrow = 1)
        
        if (pRE[g] == 1) {
          Z0 <- matrix(1, nrow = 1)
        } else if (pRE[g] == 2) {
          Z0 <- matrix(c(1, 0), nrow = 1)
        } else if (pRE[g] == 3) {
          Z0 <- matrix(c(1, 0, 0^2), nrow = 1)
        } else {
          stop("Unsupported pRE[g]")
        }
        
        suby[1, 3] <- as.numeric(Xi_0 %*% subbeta) +
          as.numeric(Z0 %*% bi[i, (index + 1):(index + pRE[g])]) +
          rnorm(1, mean = 0, sd = sqrt(sigma[g]))
        
        suby[1, 2] <- 0
        
        if (ni==0) {
          
          colnames(suby) <- c("ID", "time", "Y")
          colnames(suby)[3] <- paste0("Y", g) 
          
        } else {
          
          for (j in 1:ni) {
            
            if (pRE[g] == 1) {
              Z <- 1
            } else if (pRE[g] == 2) {
              Z <- c(1, j*increment)
            } else if (pRE[g] == 3) {
              Z <- c(1, j*increment, (j*increment)^2)
            }
            
            if (Visit[j] == 1) {
              
              time_ij <- j * increment
              Xi_t <- matrix(c(1, X[i, ], time_ij), nrow = 1)
              
              suby[j+1, 3] <- as.numeric(Xi_t %*% subbeta) +
                as.numeric(Z %*% bi[i, (index + 1):(index + pRE[g])]) +
                rnorm(1, mean = 0, sd = sqrt(sigma[g]))
              
            } else {
              suby[j+1, 3] <- NA
            }
            
            suby[j+1, 2] <- j*increment
            
          }
          
          colnames(suby) <- c("ID", "time", "Y")
          colnames(suby)[3] <- paste0("Y", g) 
          
        }
        YdataRaw <- rbind(YdataRaw, suby)
      }
      YdataRaw <- as.data.frame(YdataRaw)
      
      if(g == 1) {
        
        Ydata <- dplyr::full_join(Ydata, YdataRaw, by = "ID")
      } else {
        Ydata <- dplyr::full_join(Ydata, YdataRaw, by = c("ID", "time"))
      }
      
      index = index + pRE[g]
      
    }
    
    X <- cbind(ID, X)
    X <- as.data.frame(X)
    colnames(X)[2:3] <- c("X11", "X12")
    Ydata <- as.data.frame(Ydata)
    ydata <- dplyr::left_join(Ydata, X, by = "ID")
    cdata <- as.data.frame(cdata)
    
    a <- list(cdata, ydata)
    names(a) <- c("mvcdata", "mvydata")
    return(a)
    
  }
}
