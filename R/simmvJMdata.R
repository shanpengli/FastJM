##' Data simulation from the joint model of multivariate longitudinal biomarkers and time-to-event data
##' @title Joint modeling of multivariate longitudinal and competing risks data
##' @name simmvJMdata
##' @param seed a random seed number specified for simulating a joint model dataset.
##' @param N an integer to specify the sample size.
##' @param increment a scalar to specify the increment of visit time for longitudinal measurements.
##' @param beta a list of true parameters for the linear mixed effects sub-models. 
##' Each component of the list correspond to a specific biomarker.
##' @param sigma a vector of true error variance for all biomarkers.
##' @param gamma1 a vector of true parameters of survival fixed effects for failure 1.
##' @param gamma2 a vector of true parameters of survival fixed effects for failure 2.
##' @param alpha1 a list of true parameters for the association parameters for failure 1. 
##' Each component of the list correspond to a specific biomarker. 
##' @param alpha2 a list of true parameters for the association parameters for failure 2. 
##' Each component of the list correspond to a specific biomarker. 
##' @param lambda1 the baseline hazard rate of failure 1. 
##' An exponential distribution with a rate parameter of \code{lambda1} is assumed.
##' @param lambda2 the baseline hazard rate of failure 2. 
##' An exponential distribution with a rate parameter of \code{lambda2} is assumed.
##' @param CL a lower limit of a uniform distribution to be specified for the censoring time.
##' @param CU an upper limit of a uniform distribution to be specified for the censoring time.
##' @param covb a matrix of variance-covariance matrix of random effects.
##' @param missprob a scalar (ranging from 0 to 1) to specify the probability of 
##' missing longitudinal observations. Default is 0.
##' @param CR logical; if \code{TRUE}, simulate competing risks time-to-event data with 2 failures.
##' Default is \code{TRUE}.
##' @return  a list of datasets for both longitudinal and survival data with the elements
##' \item{mvydata}{a long-format data frame of longitudinal data.} 
##' \item{mvcdata}{a dataframe of survival data.} 
##' 
##' @export
##' 

simmvJMdata <- function(seed = 100, N = 200, increment = 0.7, beta = list(beta1 = c(5, 1.5, 2, 1),
                                                                          beta2 = c(10, 1, 2, 1)),
                        sigma = c(0.5, 0.5),
                        gamma1 = c(1, 0.5),
                        gamma2 = c(-0.5, 0.5),
                        alpha1 = list(alpha11 = c(0.5, 0.7),
                                      alpha12 = c(-0.5, 0.5)),
                        alpha2 = list(alpha21 = c(0.5, 0.7),
                                      alpha22 = c(-0.5, 0.5)),
                        lambda1 = 0.05,
                        lambda2 = 0.025,
                        CL = 5,
                        CU = 10,
                        covb = diag(rep(1, 4)),
                        missprob = 0,
                        CR = TRUE){
  
  set.seed(seed)
  
  pRE <- c()
  pREtotal = 0
  
  for(g in 1:length(beta)){
    pRE[g] <- length(alpha1[[g]])
    pREtotal = pREtotal + pRE[g]
  }
  index = 0
  
  
  if(CR == TRUE){
  bi <- MASS::mvrnorm(n = N, rep(0, length(unlist(alpha1))), covb, tol = 1e-6, empirical = FALSE)
  
  ##covariate
  X1 <- sample(c(0, 1), N, replace = TRUE, prob = c(0.5, 0.5))
  X2 <- runif(N, min = -5, max = 5)
  X <- cbind(X1, X2)
  
  #hazard rate of risk1 and risk2
  ## non-informative cencoring time 
  C <- runif(N, min = CL, max = CU)
  risk1 <- vector()
  risk2 <- vector()
  for (i in 1:N) {
    temp=lambda1*exp(X[i, ] %*% gamma1 + unlist(alpha1)%*% bi[i, ])
    risk1[i] <- rexp(1, temp)
    temp=lambda2*exp(X[i, ] %*% gamma2 + unlist(alpha2)%*% bi[i, ])
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
      
      suby[1, 3] <- subbeta[1] + subbeta[2]*X[i, 1] + subbeta[3]*X[i, 2] + 
        bi[i, index+1] + rnorm(1, mean = 0, sd = sqrt(sigma[g]))
      
      suby[1, 2] <- 0
      
      if (ni==0) {
        
        colnames(suby) <- c("ID", "time", "Y")
        colnames(suby)[3] <- paste0("Y", g) 
        
      } else {
        
        for (j in 1:ni) {
          
          if (length(alpha1[[g]]) == 2) {
            Z <- c(1, j*increment) 
          } else if (length(alpha1[[g]]) == 3) {
            Z <- c(1, j*increment, (j*increment)^2) 
          } else {
            Z <- 1
          }
          
          if (Visit[j] == 1) {
            suby[j+1, 3] <- subbeta[1] + subbeta[2]*X[i, 1] + subbeta[3]*X[i, 2] + subbeta[4]*j*increment + 
              Z%*%bi[i, (index + 1):(index + pRE[g])] + 
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
  names(a) <- c("cdata", "ydata")
  return(a)
  } else{
  
    bi <- MASS::mvrnorm(n = N, rep(0, length(unlist(alpha1))), covb, tol = 1e-6, empirical = FALSE)
  
    ##covariate
    X1 <- sample(c(0, 1), N, replace = TRUE, prob = c(0.5, 0.5))
    X2 <- runif(N, min = -5, max = 5)
    X <- cbind(X1, X2)
  
    #hazard rate of risk1 and risk2
    ## non-informative cencoring time 
    C <- runif(N, min = CL, max = CU)
    risk1 <- vector()
    for (i in 1:N) {
      temp=lambda1*exp(X[i, ] %*% gamma1 + unlist(alpha1)%*% bi[i, ])
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
    writeLines(paste0("The risk 1 rate is: ", table[2, 2], "%"))
    ID <- c(1:N)
    cdata <- cbind(ID, survtimeraw$survtime, survtimeraw$cmprsk, X)
    colnames(cdata) <- c("ID", "survtime", "cmprsk", "X21", "X22")
  
    ##fixed effects in longitudinal mean portion
    Ydata <- data.frame(c(1:N))
    colnames(Ydata) <- "ID"
    for (g in 1:length(beta)) {
    
      YdataRaw <- NULL
      subbeta <- beta[[g]]
      Z <- 1
    for (i in 1:N) {
      
      ni <- floor(cdata[i, 2]/increment)
      Visit <- sample(c(0, 1), ni, replace = TRUE, prob = c(missprob, 1 - missprob))
      suby <- matrix(0, nrow = ni+1, ncol = 3)
      suby[, 1] <- i
      
      suby[1, 3] <- subbeta[1] + subbeta[2]*X[i, 1] + subbeta[3]*X[i, 2] + 
        bi[i, index + 1] + rnorm(1, mean = 0, sd = sqrt(sigma[g]))
      
      suby[1, 2] <- 0
      
      if (ni==0) {
        
        colnames(suby) <- c("ID", "time", "Y")
        colnames(suby)[3] <- paste0("Y", g) 
        
      } else {
        
        for (j in 1:ni) {
          
          if (length(alpha1[[g]]) == 2) {
            Z <- c(1, j*increment) 
          } else if (length(alpha1[[g]]) == 3) {
            Z <- c(1, j*increment, (j*increment)^2) 
          } else {
            Z <- 1
          }
          
          if (Visit[j] == 1) {
            suby[j+1, 3] <- subbeta[1] + subbeta[2]*X[i, 1] + subbeta[3]*X[i, 2] + subbeta[4]*j*increment + 
              Z%*%bi[i,  (index + 1):(index + pRE[g])] + 
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
