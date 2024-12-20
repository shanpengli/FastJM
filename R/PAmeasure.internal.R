restricted_data_gen <- function(time, status, tau){
  output <- data.frame(time = ifelse(tau > time, time, tau),
                       status = ifelse(tau > time, status, 1))
  return(output)
}

pam.survreg <- function(fit.survreg) 
{
  x.matrix.unsorted <- fit.survreg$x
  y.unsorted <- fit.survreg$y[, 1]
  censor.unsorted <- fit.survreg$y[, 2]
  nsize <- length(y.unsorted)
  y <- sort(y.unsorted)
  delta <- censor.unsorted[order(y.unsorted)]
  p <- dim(as.matrix(x.matrix.unsorted))[2]
  if (p == 1) {
    x.matrix <- as.matrix(x.matrix.unsorted[order(y.unsorted)])
  }
  else {
    x.matrix <- x.matrix.unsorted[order(y.unsorted), ]
  }
  nsize <- length(y)
  fit.km.censoring <- survfit(Surv(y, 1 - delta) ~ 1)
  sum.km.censoring <- summary(fit.km.censoring, times = y, 
                              extend = TRUE)
  km.censoring <- sum.km.censoring$surv
  km.censoring.minus <- c(1, km.censoring[-length(km.censoring)])
  ratio.km <- delta/km.censoring.minus
  ratio.km[is.nan(ratio.km)] <- 0
  weight.km <- ratio.km/(sum(ratio.km))
  if (fit.survreg$dist == "exponential") {
    t.predicted <- predict(fit.survreg, newdata = data.frame(x.matrix), 
                           type = "response") * gamma(2)
  }
  else if (fit.survreg$dist == "weibull") {
    t.predicted <- predict(fit.survreg, newdata = data.frame(x.matrix), 
                           type = "response") * gamma(1 + fit.survreg$scale)
  }
  else if (fit.survreg$dist == "lognormal") {
    t.predicted <- predict(fit.survreg, newdata = data.frame(x.matrix), 
                           type = "response") * exp((fit.survreg$scale)^2/2)
  }
  else if (fit.survreg$dist == "loglogistic") {
    t.predicted <- predict(fit.survreg, newdata = data.frame(x.matrix), 
                           type = "response") * gamma(1 + fit.survreg$scale) * 
      gamma(1 - fit.survreg$scale)
  }
  else {
    t.predicted <- predict(fit.survreg, newdata = data.frame(x.matrix), 
                           type = "response")
  }
  wls.fitted <- tryCatch(lm(y ~ t.predicted, weights = weight.km), 
                         error = function(e) {
                           return(c(NA, NA))
                         })
  calibrate.fitted <- tryCatch(predict(wls.fitted), error = function(e) {
    return(c(NA, NA))
  })
  num.rho2 <- sum(weight.km * (calibrate.fitted - sum(weight.km * 
                                                        y))^2)
  denom.rho2 <- sum(weight.km * (y - sum(weight.km * y))^2)
  R2 <- format(round(num.rho2/denom.rho2, digits = 4), nsmall = 4)
  num.L2 <- sum(weight.km * (y - calibrate.fitted)^2)
  denom.L2 <- sum(weight.km * (y - t.predicted)^2)
  L2 <- format(round(num.L2/denom.L2, digits = 4), nsmall = 4)
  return(list(R.squared = R2, L.squared = L2))
}


predict.t.cox <- function(fit.cox) 
{
  x.matrix.unsorted <- fit.cox$x
  y.unsorted <- fit.cox$y[, 1]
  censor.unsorted <- fit.cox$y[, 2]
  my.beta <- fit.cox$coeff
  nsize <- length(y.unsorted)
  y <- sort(y.unsorted)
  delta <- censor.unsorted[order(y.unsorted)]
  p <- dim(as.matrix(x.matrix.unsorted))[2]
  if (p == 1) {
    x.matrix <- as.matrix(x.matrix.unsorted[order(y.unsorted)])
  }
  else {
    x.matrix <- x.matrix.unsorted[order(y.unsorted), ]
  }
  y.length <- length(y)
  yi.matrix <- matrix(rep(y, each = y.length), nrow = y.length)
  yj.matrix <- t(yi.matrix)
  fit.km.censoring <- survfit(Surv(y, 1 - delta) ~ 1)
  sum.km.censoring <- summary(fit.km.censoring, times = y, 
                              extend = TRUE)
  km.censoring <- sum.km.censoring$surv
  km.censoring.minus <- c(1, km.censoring[-length(km.censoring)])
  ratio.km <- delta/km.censoring.minus
  ratio.km[is.nan(ratio.km)] <- 0
  weight.km <- ratio.km/(sum(ratio.km))
  R <- ((yj.matrix >= yi.matrix) * 1)
  my.Lambda <- R %*% ((delta)/t(t(exp(x.matrix %*% my.beta)) %*% 
                                  R))
  rm(R)
  my.power <- matrix(rep(t(exp(x.matrix %*% my.beta)), each = y.length), 
                     nrow = y.length)
  my.factor <- (max(my.Lambda)/100)
  my.Lambda2 <- t(matrix(rep(exp(-my.Lambda/my.factor), each = y.length), 
                         nrow = y.length))
  S.hat.x <- my.Lambda2^(my.factor * my.power)
  rm(my.Lambda)
  rm(my.power)
  rm(my.Lambda2)
  t1 <- y
  t2 <- c(0, t1[1:length(t1) - 1])
  delta.t <- t1 - t2
  t.predicted <- colSums(delta.t * S.hat.x)
  
  return(data.frame(y = y, t.predicted = t.predicted, 
                    risk.score = x.matrix %*% my.beta,
                    status = ifelse(delta==1, "event", "censored"),
                    assumption = "cox"))
}

restricted_pa_cox <- function(fit.cox, y.input = NULL, predict = F) 
{
  x.matrix.unsorted <- fit.cox$x
  y.unsorted <- fit.cox$y[, 1]
  y.order <- order(y.unsorted)
  censor.unsorted <- fit.cox$y[, 2]
  my.beta <- fit.cox$coeff
  nsize <- length(y.unsorted)
  y <- y.unsorted[y.order]
  delta <- censor.unsorted[y.order]
  p <- dim(as.matrix(x.matrix.unsorted))[2]
  if (p == 1) {
    x.matrix <- as.matrix(x.matrix.unsorted[y.order])
  }
  else {
    x.matrix <- x.matrix.unsorted[y.order, ]
  }
  y.length <- length(y)
  yi.matrix <- matrix(rep(y, each = y.length), nrow = y.length)
  yj.matrix <- t(yi.matrix)
  fit.km.censoring <- survfit(Surv(y, 1 - delta) ~ 1)
  sum.km.censoring <- summary(fit.km.censoring, times = y, 
                              extend = TRUE)
  km.censoring <- sum.km.censoring$surv
  km.censoring.minus <- c(1, km.censoring[-length(km.censoring)])
  ratio.km <- delta/km.censoring.minus
  ratio.km[is.nan(ratio.km)] <- 0
  weight.km <- ratio.km/(sum(ratio.km))
  R <- ((yj.matrix >= yi.matrix) * 1)
  my.Lambda <- R %*% ((delta)/t(t(exp(x.matrix %*% my.beta)) %*% 
                                  R))
  rm(R)
  my.power <- matrix(rep(t(exp(x.matrix %*% my.beta)), each = y.length), 
                     nrow = y.length)
  my.factor <- (max(my.Lambda)/100)
  my.Lambda2 <- t(matrix(rep(exp(-my.Lambda/my.factor), each = y.length), 
                         nrow = y.length))
  S.hat.x <- my.Lambda2^(my.factor * my.power)
  rm(my.Lambda)
  rm(my.power)
  rm(my.Lambda2)
  
  t1 <- y.input[y.order]
  t2 <- c(0, t1[1:length(t1) - 1])
  delta.t <- t1 - t2
  t.predicted <- colSums(delta.t * S.hat.x)
  
  if(predict == TRUE){return(data.frame(pred = t.predicted, 
                                        obs = t1,
                                        status = delta))}
  
  y <- y.input[y.order]
  
  wls.fitted <- tryCatch(lm(y ~ t.predicted, weights = weight.km), 
                         error = function(e) {
                           return(c(NA, NA))
                         })
  calibrate.fitted <- tryCatch(predict(wls.fitted), error = function(e) {
    return(c(NA, NA))
  })
  num.rho2 <- sum(weight.km * (calibrate.fitted - sum(weight.km * 
                                                        y))^2)
  denom.rho2 <- sum(weight.km * (y - sum(weight.km * y))^2)
  R2 <- round(num.rho2/denom.rho2, digits = 4)
  num.L2 <- sum(weight.km * (y - calibrate.fitted)^2)
  denom.L2 <- sum(weight.km * (y - t.predicted)^2)
  L2 <- round(num.L2/denom.L2, digits = 4)
  SR <- round(R2 * L2, digits = 4)
  return(list(R.squared = format(R2, nsmall = 4), 
              L.squared = format(L2, nsmall = 4), 
              Psuedo.R = format(SR, nsmall = 4)))
}


# aft models
restricted_pa_aft <- function(fit.survreg, tau, y.input = NULL, predict = F) 
{
  x.matrix.unsorted <- fit.survreg$x
  y.unsorted <- fit.survreg$y[, 1]
  censor.unsorted <- fit.survreg$y[, 2]
  nsize <- length(y.unsorted)
  y.order <- order(y.unsorted)
  y <- y.unsorted[y.order]
  delta <- censor.unsorted[y.order]
  p <- dim(as.matrix(x.matrix.unsorted))[2]
  if (p == 1) {
    x.matrix <- as.matrix(x.matrix.unsorted[y.order])
  }
  else {
    x.matrix <- x.matrix.unsorted[y.order, ]
  }
  nsize <- length(y)
  fit.km.censoring <- survfit(Surv(y, 1 - delta) ~ 1)
  sum.km.censoring <- summary(fit.km.censoring, times = y, 
                              extend = TRUE)
  km.censoring <- sum.km.censoring$surv
  km.censoring.minus <- c(1, km.censoring[-length(km.censoring)])
  ratio.km <- delta/km.censoring.minus
  ratio.km[is.nan(ratio.km)] <- 0
  weight.km <- ratio.km/(sum(ratio.km))
  
  if (fit.survreg$dist %in% c("exponential", "weibull")) {
    temp.pred <- predict(fit.survreg, newdata = data.frame(x.matrix), 
                         type = "response")
    gtau <- (tau/temp.pred)^(1/fit.survreg$scale)
    t.predicted <- temp.pred * gamma(1 + fit.survreg$scale) *
      (1-gammainc((1 + fit.survreg$scale), gtau))
  }
  else if (fit.survreg$dist == "lognormal") {
    temp.pred <- predict(fit.survreg, newdata = data.frame(x.matrix), 
                         type = "response")
    gtau <- (log(tau)-log(temp.pred))/fit.survreg$scale
    t.predicted <- temp.pred * exp((fit.survreg$scale)^2/2) * pnorm(gtau)
  }
  
  y <- y.input[y.order]
  
  if(predict == T){return(data.frame(pred = t.predicted, 
                                     obs = y,
                                     status = delta))}
  
  wls.fitted <- tryCatch(lm(y ~ t.predicted, weights = weight.km), 
                         error = function(e) {
                           return(c(NA, NA))
                         })
  calibrate.fitted <- tryCatch(predict(wls.fitted), error = function(e) {
    return(c(NA, NA))
  })
  num.rho2 <- sum(weight.km * (calibrate.fitted - sum(weight.km * 
                                                        y))^2)
  denom.rho2 <- sum(weight.km * (y - sum(weight.km * y))^2)
  R2 <- round(num.rho2/denom.rho2, digits = 4)
  num.L2 <- sum(weight.km * (y - calibrate.fitted)^2)
  denom.L2 <- sum(weight.km * (y - t.predicted)^2)
  L2 <- round(num.L2/denom.L2, digits = 4)
  SR <- round(R2 * L2, digits = 4)
  return(list(R.squared = format(R2, nsmall = 4), 
              L.squared = format(L2, nsmall = 4), 
              Psuedo.R = format(SR, nsmall = 4)))
}

compute_cindex <- function(predicted_time, survival_time,
                           status, weight = "H", 
                           input_tau = NULL) {
  ties <- 0
  concordant <- 0
  usable_pairs <- 0
  if(weight == "U"){
    KM <- summary(survfit(Surv(survival_time, status==0)~1))
    df_KM <- data.frame(time = c(0, KM$time), surv = c(1, KM$surv))
    for (i in 1:length(survival_time)) {
      for (j in 1:length(survival_time)) {
        if(status[i] == 1){
          usable_pairs <- usable_pairs + 
            (survival_time[i] < survival_time[j]) * 
            min(df_KM$surv[df_KM$time <= survival_time[i]])^(-2)
          concordant <- concordant +
            (predicted_time[i] < predicted_time[j]) * 
            (survival_time[i] < survival_time[j]) * 
            min(df_KM$surv[df_KM$time <= survival_time[i]])^(-2)
        }
      }
    }
  } else if(weight == "H"){
    for (i in 1:length(survival_time)) {
      for (j in 1:length(survival_time)) {
        if(i != j & status[i] == 1){
          ties <- ties + 
            (predicted_time[i] == predicted_time[j]) *
            (survival_time[i] != survival_time[j])
          usable_pairs <- usable_pairs + 
            (survival_time[i] < survival_time[j]) * 
            (predicted_time[i] != predicted_time[j])
          concordant <- concordant +
            (predicted_time[i] < predicted_time[j]) * 
            (survival_time[i] < survival_time[j])
        }
      }
    }
  } else if(weight == "U_tau"){
    KM <- summary(survfit(Surv(survival_time, status==0)~1))
    df_KM <- data.frame(time = c(0, KM$time), surv = c(1, KM$surv))
    for (i in 1:length(survival_time)) {
      for (j in 1:length(survival_time)) {
        if(status[i] == 1){
          ties <- ties + 
            (predicted_time[i] == predicted_time[j]) *
            (survival_time[i] != survival_time[j])*
            (survival_time[i] < input_tau)*
            min(df_KM$surv[df_KM$time <= survival_time[i]])^(-2)*
            (i <= j)
          usable_pairs <- usable_pairs + 
            (survival_time[i] < survival_time[j]) * 
            (survival_time[i] < input_tau) *
            min(df_KM$surv[df_KM$time <= survival_time[i]])^(-2)* 
            (predicted_time[i] != predicted_time[j])
          concordant <- concordant +
            (predicted_time[i] < predicted_time[j]) * 
            (survival_time[i] < survival_time[j]) *
            (survival_time[i] < input_tau) *
            min(df_KM$surv[df_KM$time <= survival_time[i]])^(-2)
        }
      }
    }
  } 
  print(ties)
  print(concordant)
  print(usable_pairs)
  cindex <- (concordant + ties/2) / (usable_pairs + ties)
  #print(c(concordant, usable_pairs))
  return(cindex)
}


##### R^2_{sph}
##### sourced from: 
##### https://ibmi3.mf.uni-lj.si/ibmi-english/biostat-center/programje/Re.r
#includes the following functions:

#re ... function for calculating the Re measure
#called with
#re(fit)
# works with coxph, survreg and aareg

#summary.re ... function for calculating the Re in time (cumulative and time window)
#called with
#summary(re(fit))


re <- function (fit, ...) UseMethod("re")



re.coxph <- function(fit,Gmat){
  #################################################################################################
  # A function that calculates the Re measure of explained variation of the model
  # Arguments: fit: A Cox model object
  #	     Gmat: a matrix of weights - if omitted, calculated as inverse KM
  # 
  # Returned values:
  #  Re: the Re measure
  #  Re.imp: the Re measure assuming constant coefficients after the last event time
  #  Re.fix: the Re measure asuming a constant value of coefficient throughout the follow-up
  #          time(gives NA in case of counting process type times (start,stop)) 
  #      se: standard error of re
  #    sen0: standard error of re, calculated with expected Q	
  #       C: C index (generalized to allow for time-dependent covariates and weighted to ensure independence of 
  #	   the censoring process)
  #    r2nw: unweighted measure
  #
  # Notes: left truncated data not yet implemented.
  #	 when calculating the re.imp, the data-set should not be split after the last event time, i.e. each
  #	 data line should represent a different subject
  #	
  #################################################################################################
  
  Y <- fit$y
  f.type <- attr(Y, "type")  
  if (f.type== "right") Y <- cbind(rep(0,nrow(Y)),Y)
  sort.it <- order(Y[,2],-Y[,3]) 
  bx <- fit$linear.predictors[sort.it] - mean(fit$linear.predictors)
  Y <- Y[sort.it,]
  
  n.ind <- nrow(Y)
  ti <- sort(unique(Y[Y[,3]==1,2]))			#ordered unique times of events
  n.times <- length(ti)	 				#number of unique times of events	
  
  
  r2<-NULL
  if(missing(Gmat)){
    km.inv <- my.survfit(Y[,1],Y[,2],Y[,3])		#now includes time-dependent data
    G <- km.inv$surv.i2[km.inv$n.event != 0]								
    dSt <- 1/G						#is this needed (we already have only at event times)
    Gmat <- matrix(G,byrow=T,ncol=length(G),nrow=n.ind)	#times in columns, ind. in rows
  }
  else if(is.null(dim(Gmat))){
    if(length(Gmat)!=n.times)stop("Wrong length of weights")
    else Gmat <- matrix(Gmat,byrow=T,ncol=length(Gmat),nrow=n.ind)
  }
  else if(any(dim(Gmat)!=c(n.ind,n.times)))stop("Wrong dimension of weights")
  Gmat <- Gmat[sort.it,]					#sort with respect to event times
  
  
  ebxun <- ebx <- exp(bx)
  rg <- meanr <- smin <- newvar1 <- newvar2 <- newvar3 <- eden <- enum <-  rep(NA,n.times)
  for(it in 1:length(ti)){
    k <- sum(Y[Y[,3]==1,2]==ti[it])			#number of ties
    
    inx <- (Y[,2] >= ti[it] & Y[,1]< ti[it])	#at risk
    #reverse the order to get the right ranks
    rangi <- (rank(-bx[inx]) - .5)/Gmat[inx,it]+.5	#ranks for each individual
    
    rg[it] <- sum(rangi[1:k]/(Gmat[inx,it])[1:k])	#rank of the one who had the event (sum of ranks if ties)
    
    r0 <- mean(rangi)
    meanr[it] <-  r0*sum(1/((Gmat[inx,it])[1:k]))	#mean rank at this time 
    #browser()
    rP <- (.5+ .5*sum((Y[inx,2] == ti[it] & Y[inx,1]< ti[it])/Gmat[inx,it]))
    smin[it] <- rP*sum(1/((Gmat[inx,it])[1:k]))	#min rank (1 if no ties)
    
    pi <- ebx[inx]/sum(ebx[inx])			#individual hazard (breslow * ebx)
    newvar1[it] <- sum((rangi-r0)^2*pi/(Gmat[inx,it])^2)	#term1 -  for nom
    newvar2[it] <- sum((rP-r0)^2*pi/(Gmat[inx,it])^2)	#term2 -  for denom
    newvar3[it] <- sum((r0-rangi)*(r0-rP)*pi/(Gmat[inx,it])^2)	#term3 - covariance
    enum[it] <- sum((r0-rangi)*pi/(Gmat[inx,it]))		#expected numerator=Q1
    eden[it] <- sum((r0-rP)*pi/(Gmat[inx,it]))	#expected denominator=Q2
    
  }
  num <- sum(meanr-rg)
  den<- sum(meanr-smin)
  eden <- sum(eden)
  enum <- sum(enum)
  r2 <- num/den	#the Re measure
  r2nw <- sum((meanr-rg))/sum((meanr-smin))
  newvar <- sqrt(sum(newvar1)/den^2-2*sum(newvar3)*num/den^3+sum(newvar2)*num^2/den^4)
  newvare <- sqrt(sum(newvar1)/eden^2-2*sum(newvar3)*enum/eden^3+sum(newvar2)*enum^2/eden^4)
  
  ebx <- sort(ebx)					#this is for Jfixed
  Jf <- NA
  if(f.type=="right"){
    nd <- length(ebx)
    imat <- matrix(ebx,nrow=nd,ncol=nd)
    jmat <- matrix(ebx,nrow=nd,ncol=nd,byrow=T)
    iimat <- imat/(imat+jmat)
    jjmat <- jmat/(imat+jmat)
    Jf <- 2*sum((jjmat-iimat)[upper.tri(jjmat)])/(sum(upper.tri(jjmat))*2)
  }
  
  #and now for the part after last event time:
  ebxun <- ebxun[Y[,2]>=max(ti)&Y[,3]==0]
  
  km.inv <- my.survfit(Y[,1],Y[,2],Y[,3])		
  nd <- sum(km.inv$n.cens[km.inv$time>max(ti)])
  if((nd>0)&(!any(km.inv$n.cens[km.inv$time>max(ti)]==0))){		#if there is at least 1 indiv. and no time-splitting
    
    ebx <- sort(ebxun)
    sumt <- 0
    for(it in 1:nd){
      sumt1 <- sum((ebx[it]-ebx[1:it])/(ebx[it]+ebx[1:it]))
      sumt <- sumt + sumt1
    }
    Jfin <- sumt/2
    
    Jfin <- (num+Jfin)/(den +  (nd*(nd-1))/4)
  }
  
  else Jfin <- NA 		
  
  out <- list(Re=r2,Re.imp=Jfin,Re.fix=Jf,se=newvare,C=(r2+1)/2,sen0=newvar,times=ti,meanr=meanr,ranks=rg,weights=Gmat,perfr=smin,type="coxph",nd=nd,r2nw=r2nw)
  class(out) <- "re"
  out
}

re.aareg <- function(fit,Gmat){
  #################################################################################################
  # A function that calculates the Re measure of explained variation of the model
  # Arguments: fit: An aareg model object
  #	     Gmat: a matrix of weights - if omitted, calculated as inverse KM
  # 
  # Returned values:
  #  Re: the Re measure
  #  
  #      se: standard error of re
  #       C: C index (generalized to allow for time-dependent covariates and weighted to ensure independence of 
  #	   the censoring process)
  # 
  #
  # Notes: left truncated data not yet implemented.
  #	
  #################################################################################################
  
  if(length(fit$y)==0)stop("The aareg model must be fitted with 'y=TRUE' option")
  if(length(fit$x)==0)stop("The aareg model must be fitted with 'x=TRUE' option")
  Y <- fit$y
  sort.it <- order(Y[,2],-Y[,3]) 
  X <- fit$x
  X <- cbind(rep(1,nrow(X)),X)
  Y[Y[,1]==-1,1] <- 0
  Y <- Y[sort.it,]
  X <- X[sort.it,,drop=FALSE]
  n.ind <- nrow(Y)
  
  ti <- sort(unique(Y[Y[,3]==1,2]))			#ordered unique times of events
  n.times <- length(ti)	 				#number of unique times of events	
  
  #nelson-Aalen:
  d <- survfit(Surv(Y[,1],Y[,2],Y[,3])~1,type="fleming-harrington")
  cumh <- -log(d$surv[d$n.event!=0])
  haz <- c(cumh[1],diff(cumh))
  
  #change the output for ties in surv. time (sum up the contributions in those intervals)
  d1 <- fit$times
  d2 <- c(d1[-1],d1[length(d1)]+1)
  ic <- which(d1!=d2)
  fc <- apply(fit$coef,2,cumsum)
  fc <- fc[ic,]
  fd <- apply(fc,2,function(x)c(x[1],diff(x)))
  #add lines for times at the end	
  
  if(n.times>nrow(fd)){
    nman <- n.times-nrow(fd)
    fd <- rbind(fd,matrix(0,nrow=nman,ncol=ncol(fd)))		
    fd[(n.times-nman+1):n.times,1] <- haz[(n.times-nman+1):n.times]
  }
  
  bx <- (X%*%t(fd))
  
  if(missing(Gmat)){
    km.inv <- my.survfit(Y[,1],Y[,2],Y[,3])		#now includes time-dependent data
    G <- km.inv$surv.i2[km.inv$n.event != 0]							
    dSt <- 1/G
    Gmat <- matrix(G,byrow=T,ncol=length(G),nrow=n.ind)	#casi v stolpcih, ljudje v vrsticah	
  }
  else if(is.null(dim(Gmat))){
    if(length(Gmat)!=n.times)stop("Wrong length of weights")
    else Gmat <- matrix(Gmat,byrow=T,ncol=length(Gmat),nrow=n.ind)
  }
  else if(any(dim(Gmat)!=c(n.ind,n.times)))stop("Wrong dimension of weights")
  
  
  rg <- meanr <- smin <- newvar1 <- newvar2 <- newvar3 <- enum <- eden <-  rep(NA,n.times)
  for(it in 1:n.times){
    
    k <- sum(Y[Y[,3]==1,2]==ti[it])			#number of ties
    
    inx <- (Y[,2] >= ti[it] & Y[,1]< ti[it])	#at risk
    pi <- bx[inx,it]			#hazard for each indiv.
    
    rangi <- (rank(-bx[inx,it]) - .5)/Gmat[inx,it]+.5	#ranks for each individual
    
    r0 <- mean(rangi)
    rP <- (.5+ .5*sum((Y[inx,2] == ti[it] & Y[inx,1]< ti[it])/Gmat[inx,it]))
    
    rg[it] <- sum(rangi[1:k]/(Gmat[inx,it])[1:k])	#rank of the one who had the event (sum of ranks if ties)
    meanr[it] <-  r0*sum(1/((Gmat[inx,it])[1:k]))	#mean rank at this time 
    smin[it] <- rP*sum(1/((Gmat[inx,it])[1:k]))		
    
    newvar1[it] <- sum((rangi-r0)^2*pi/(Gmat[inx,it])^2)	#term1 -  for nom
    newvar2[it] <- sum((rP-r0)^2*pi/(Gmat[inx,it])^2)	#term2 -  for denom
    newvar3[it] <- sum((r0-rangi)*(r0-rP)*pi/(Gmat[inx,it])^2)	#term3 - covariance
    enum[it] <- sum((r0-rangi)*pi/(Gmat[inx,it]))		#expected numerator
    eden[it] <- sum((r0-rP)*pi/(Gmat[inx,it]))	#expected denominator
    
  }
  num <- sum(meanr-rg)
  den<- sum(meanr-smin)
  enum <- sum(enum)
  eden <- sum(eden)
  r2 <- num/den	#the Re measure
  newvar <- sqrt(sum(newvar1)/den^2-2*sum(newvar3)*num/den^3+sum(newvar2)*num^2/den^4)
  newvare <- sqrt(sum(newvar1)/eden^2-2*sum(newvar3)*enum/eden^3+sum(newvar2)*enum^2/eden^4)
  out <- list(Re=r2,se=newvare,C=(r2+1)/2,sen0=newvar,type="aareg")
  class(out) <- "re"
  out
}

re.survreg <- function(fit,Gmat){
  #################################################################################################
  # A function that calculates the Re measure of explained variation of the model
  # Arguments: fit: A survreg model object
  #	    
  # 
  # Returned values:
  #  Re: the Re measure
  #      se: standard error of re
  #       C: C index (generalized to allow for time-dependent covariates and weighted to ensure independence of 
  #	   the censoring process)
  #  r2nw,senw: unweighted measure
  #
  #Notes: times should not be negative
  #	
  #################################################################################################
  
  Y <- fit$y
  Y <- cbind(rep(0,nrow(Y)),Y)
  if(any(Y[,2]<0))warning("Negative follow-up times")
  #nazaj antilogaritmiramo cas:
  #if(fit$dist%in%c("weibull","exponential","loglogistic","lognormal"))Y[,2] <- exp(Y[,2])
  
  sort.it <- order(Y[,2],-Y[,3]) 
  Y <- Y[sort.it,]
  
  if(fit$dist=="weibull"|fit$dist=="exponential")bx <- -(fit$linear.predictors[sort.it] - fit$coef[1])
  bx <- -fit$linear.predictors[sort.it]
  if(fit$dist=="weibull"|fit$dist=="exponential"){
    gamma <- (1/fit$scale)
    lambda <- exp(-fit$coef[1])^gamma
  }
  
  n.ind <- nrow(Y)
  ti <- sort(unique(Y[Y[,3]==1,2]))			#ordered unique times of events
  n.times <- length(ti)	 				#number of unique times of events	
  
  r2<-NULL
  
  if(missing(Gmat)){
    km.inv <- my.survfit(Y[,1],Y[,2],Y[,3])		#now includes time-dependent data
    G <- km.inv$surv.i2[km.inv$n.event != 0]							
    dSt <- 1/G
    Gmat <- matrix(G,byrow=T,ncol=length(G),nrow=n.ind)	#casi v stolpcih, ljudje v vrsticah	
  } else if(is.null(dim(Gmat))){
    if(length(Gmat)!=n.times)stop("Wrong length of weights")
    else Gmat <- matrix(Gmat,byrow=T,ncol=length(Gmat),nrow=n.ind)
  } else if(any(dim(Gmat)!=c(n.ind,n.times)))stop("Wrong dimension of weights")
  
  
  rg <- meanr <- smin <- newvar1 <- newvar2 <- newvar3 <- eden <- enum  <- rep(NA,n.times)
  for(it in 1:length(ti)){
    k <- sum(Y[Y[,3]==1,2]==ti[it])			#number of ties
    inx <- (Y[,2] >= ti[it] & Y[,1]< ti[it])	#at risk
    if(fit$dist=="weibull"|fit$dist=="exponential"){
      if(it==1) pi <- lambda*exp(bx)[inx]*(ti[it])^gamma						
      else pi <- lambda*exp(bx)[inx]*((ti[it])^gamma-(ti[it-1])^gamma)	
    }
    else if(fit$dist=="gaussian"|fit$dist=="lognormal"){
      bxi <- -bx[inx]		#dam ga nazaj na pozitivne
      casi <- ti
      if(fit$dist=="lognormal") {
        casi <- log(ti)
      }
      S  <- 1-pnorm((casi[it]- bxi)/fit$scale)
      S0 <- 1
      if(it!=1) S0 <- 1-pnorm((casi[it-1]- bxi)/fit$scale)
      pi <- -log(S/S0)
    }			
    else if(fit$dist=="logistic"|fit$dist=="loglogistic"){
      bxi <- -bx[inx]
      casi <- ti
      if(fit$dist=="loglogistic") {
        casi <- log(ti)
      }
      w <- exp((casi[it]- bxi)/fit$scale)
      S <- 1/(1+w)
      S0 <- 1
      if(it!=1){
        w0 <- exp((casi[it-1]- bxi)/fit$scale)
        S0 <- 1/(1+w0)
      }
      pi <- -log(S/S0)
    }			
    rangi <- (rank(-bx[inx]) - .5)/Gmat[inx,it]+.5	#ranks for each individual
    
    r0 <- mean(rangi)
    rP <- (.5+ .5*sum((Y[inx,2] == ti[it] & Y[inx,1]< ti[it])/Gmat[inx,it]))
    
    rg[it] <- sum(rangi[1:k]/(Gmat[inx,it])[1:k])	#rank of the one who had the event (sum of ranks if ties)
    meanr[it] <-  r0*sum(1/((Gmat[inx,it])[1:k]))	#mean rank at this time 
    smin[it] <- rP*sum(1/((Gmat[inx,it])[1:k]))		
    
    newvar1[it] <- sum((rangi-r0)^2*pi/(Gmat[inx,it])^2)	#term1 -  for nom
    newvar2[it] <- sum((rP-r0)^2*pi/(Gmat[inx,it])^2)	#term2 -  for denom
    newvar3[it] <- sum((r0-rangi)*(r0-rP)*pi/(Gmat[inx,it])^2)	#term3 - covariance
    enum[it] <- sum((r0-rangi)*pi/(Gmat[inx,it]))		#expected numerator
    eden[it] <- sum((r0-rP)*pi/(Gmat[inx,it]))	#expected denominator
    
  }
  num <- sum(meanr-rg)
  den<- sum(meanr-smin)
  enum <- sum(enum)
  eden <- sum(eden)
  r2 <- num/den	#the Re measure
  newvar <- sqrt(sum(newvar1)/den^2-2*sum(newvar3)*num/den^3+sum(newvar2)*num^2/den^4)
  newvare <- sqrt(sum(newvar1)/eden^2-2*sum(newvar3)*enum/eden^3+sum(newvar2)*enum^2/eden^4)
  
  out <- list(Re=r2,se=newvare,C=(r2+1)/2,sen0=newvar,times=ti,meanr=meanr,perfr=smin,ranks=rg,weights=Gmat,type="survreg")
  class(out) <- "re"
  out	
}



print.re <- function(x, digits=4, ...){
  
  cat("Re measure=",round(x$Re,digits), ",    SE=", round(x$se,digits), "\n \n",sep="")
  if(x$type=="coxph"){
    if(x$nd){
      cat("Number of censored after the last failure time: ", x$nd, "\n")
      if(!is.na(x$Re.imp))cat("Re corrected for the censoring at the end \n (assuming the coeff. and covariate values fixed after the last failure time): ", round(x$Re.imp,digits),"\n \n")
      else cat("Re corrected for the censoring at the end can not be computed - the follow-up \n time after the last event time is split into several intervals \n")
    }
    if(!is.na(x$Re.fix)){
      cat("Integrated version of Re (assuming all coeff. and covariate values fixed):", round(x$Re.fix,digits),"\n")	
    }
  }
  cat("\n")
  cat("Generalized C index=",round(x$C,digits), "\n")
  invisible(x)
}



summary.re <- function(object, times, band=5,...){
  
  ####################################
  # A function for calculating the RE in time 
  # object ... an result of a re(fit) call
  # times ... times at which we wish to calculate Re (unique event times by default
  # band ... a window for calculating the dRe 
  #######################################
  
  num <- cumsum(object$meanr-object$ranks)
  den <- cumsum(object$meanr-object$perfr) 
  
  Rti <- num/den
  
  dn <- length(num) - band 
  
  num1 <- num[-(1:band)]
  den1 <- den[-(1:band)]
  
  dRti <- (num1 - num[1:dn])/(den1 - den[1:dn])
  
  dRti <- c(rep(NA,band-1),Rti[band],dRti)
  
  tab <- data.frame(times=object$times,Rti=Rti,dRti=dRti)
  if(!missing(times)) {
    search.time <- function(time,times)max(which(times < time))
    inx <- unlist(lapply(times,search.time, tab$times)	)
    tab[inx,1] <- times
    tab <- tab[inx,]
  }
  tab
}
























"my.survfit" <- function(start,stop,event){
  #internal function for calculation of weights (allows counting process notation)
  
  if(missing(start)) start <- rep(0,length(stop))		#if no count. process notation, set everyone to start at 0
  
  data <- data.frame(start=start,stop=stop,event=event)	
  
  data <- data[order(data$start),]			#ordered according to the starting time
  data$Y <- unlist(lapply(data$stop,function(x)sum(data$start<=x)))	#for each time point: the number of indiv. that started before (or at) it
  
  data <- data[order(data$stop,-data$event),]		#order in time, put events before censoring
  
  howmany <- nrow(data)					#total number of rows
  
  lt.t1 <- c(0,data$stop[-nrow(data)])			#aux. vector for searching the stop times	
  lt.t2 <- c(data$stop[-1],data$stop[nrow(data)]+1)
  lt1 <- which(data$stop!=lt.t1)				#first case at a time point
  lt2 <- which(data$stop!=lt.t2)				#last case at a time point
  
  ct.t1 <- c(-1,data$start[-nrow(data)])			#aux. vector for searching the start times
  ct.t2 <- c(data$start[-1],data$start[nrow(data)]+1)
  ct1 <- which(data$start!=ct.t1)				#first case at a time point
  ct2 <- which(data$start!=ct.t2)				#last case at a time point
  n.incom <- rep(0,length(lt1))				#prepare a vector that will contain the no. of incoming indiv. at a time point (length=no. of unique time points)
  inx <- unlist(lapply(data$start[ct1],function(x)min(which(data$stop[lt1]>=x)))) #the index of the first one to fail after a starting point
  incom.sum <- diff(c(0,ct2))				#the number of incoming
  if(length(inx)>1){
    if(inx[1]==inx[2]){					#if indiv. come in at time 0 and first stop time
      incom.sum <- c(sum(incom.sum[1:2]),incom.sum[-(1:2)])    #sum them up
      inx <- unique(inx)					#only one value per each stop time
    }
  }
  n.incom[inx] <- incom.sum				#at each event time, we need the number of patients that got in just before it.
  n.incom2 <- n.incom					
  n.incom2[1] <- n.incom2[1] - sum(data$start==0)		#the number of additional patients to come in at time 0 is set to 0
  
  
  evcum <- cumsum(data$event)			#number of events up to a certain time point
  evcum <- evcum[lt2]				#number of events for t<=t_i
  cencum <- cumsum(data$event==0)			#number of censored
  cencum <- cencum[lt2]				#number of censored t<=t_i
  n.event <- diff(c(0,evcum))			#number of events at each unique event time
  n.cens <- diff(c(-sum(data$start==0),cencum))	#number of censored at each unique event time, the first number is the total of the patients
  n.t.cens <- n.cens - n.incom			#number truly censored: number censored - number of incoming (0 at time 0)
  n.risk <- (data$Y - (0:(howmany-1)))[lt1]-n.incom2 #number at risk at each event time: no.started before or at this time - no. of incoming at this time - 1 per row before this time
  
  kmji <- 1- n.event/n.risk			#proportion of survivors
  
  kmji.i <- 1- n.t.cens/(n.risk-n.event)		#the events are taken out of the risk set before calculating the G(t_i)
  
  km <- (cumprod(kmji))				#km	
  km.i <- cumprod(kmji.i)				#km for censoring			
  
  #correction for the last value (if division by 0) (this happens if the last time is an event time)
  inx <- which(is.na(km.i))				#missing at the end, only needed if not lagged?
  if(length(inx)){
    km.i[inx] <- km.i[min(inx)-1]	# carrie last value forward
  }
  
  km.i2 <- c(1,km.i[-length(km.i)])			#this gives the lagged weight (weight at time just before t)
  
  list(surv=km,surv.i=km.i,n.event=n.event,n.cens=n.t.cens,time=data$stop[lt1],n.risk=n.risk,surv.i2=km.i2)
}

############## 2024-8-7 add competing risk model

##################################################
## calculate integral for int_0^tau t dF(t|x) 
##################################################
### input:
# x: input estimated CIF for subject i
# time.cif: time point of observation event 1
# tau: restricted time
m_cif <- function(x, time.cif, tau){
  index <- which(time.cif <= tau)
  t(diff(c(0, time.cif[index], tau))) %*% c(1, (1-x)[index])
  #t(c(time.cif[index], tau))%*%(diff(c(0, x[index], 1)))
}

##################################################
### calculate PA measure for competing risk model
##################################################
### input:
# ftime: vector of failure/censoring times
# fstatus: vector with a unique code for each failure type (1,2,3...)
#          and censored observations (0)
# tau: restricted time
# time.cif: time point of observation event of interest
# pred.cif: estimated CIF for all subject
pam.censor.cr <- function(ftime, fstatus, tau, 
                          pred.cif, time.cif,
                          event.type = 1){
  
  # create new variable T(1,tau)
  ftime.new <- ifelse(ftime > tau, tau, ftime)
  ftime.new <- ifelse(ftime.new <= tau & !(fstatus %in% c(event.type, 0)), 
                      tau, ftime.new)
  fstatus.new <- ifelse(ftime.new == tau | fstatus > 0, 1, fstatus)
  
  # estimate m(1,tau) 
  pred <- apply(pred.cif, 2, m_cif, time.cif = time.cif, tau = tau)
  
  # obtain PA measure
  out <- pam.censor(ftime.new, pred, fstatus.new)
  out <- c(out, 
           Psuedo.R = format(round(as.numeric(out$R.squared) * 
                                     as.numeric(out$L.squared), 
                                   digits = 4), nsmall = 4))
  return(out)
}

######################################################
### calculate baseline cause specific hazard rates h0
######################################################
### input:
# ftime: vector of failure/censoring times
# fstatus: vector with a unique code for each failure type (1,2,3...)
#          and censored observations (0)
# cox.model: cox model applied on the target cause
### output: 
# baseline cause specific hazard rates h0
# sorted event time y
get_h0 <- function(ftime, fstatus, cox.model, X){
  tab <- data.frame(table(ftime[fstatus == 1])) 
  y <- as.numeric(levels(tab[, 1]))[tab[, 1]]
  d <- tab[, 2]
  betaHat <- cox.model$coef
  h0 <- rep(NA, length(y))
  for(l in 1:length(y)){
    h0[l] <- d[l] / sum(exp(X[ftime >= y[l], ] %*% betaHat))
  }
  return(data.frame(h0, y))
}

####################################################################
### calculate CIF based on baseline cause specific hazard rates h0
####################################################################
### input:
# ftime: vector of failure/censoring times
# fstatus: vector with a unique code for each failure type (1,2,3...)
#          and censored observations (0)
# cox.focus: cox model applied on the target cause
# cox.other: cox models applied on other cause (store all other cause model in a list)
# X: Xariates input for models
### output: 
# sorted event time y
# CIF

get_CIF <- function(ftime, fstatus, cox.focus, cox.other, X){
  
  fit <- tryCatch({survfit(cox.focus, newdata = X, se.fit = FALSE)}, 
                  error = handle_error)
  Ht <- fit$cumhaz[fit$n.event==1,]
  if(is.null(dim(Ht))) return()
  hi <- apply(Ht, 2, function(x) diff(c(0, x)))
  
  for(cox.model in cox.other){
    fit_other <- tryCatch({survfit(cox.model, newdata = X, se.fit = FALSE)}, 
                          error = handle_error)
    Ht_other <- fit_other$cumhaz[fit$n.event==1,]
    Ht <- Ht + Ht_other
    if(is.null(dim(Ht))) return()
  }
  
  CIF <- apply(exp(-Ht) * hi, 2, cumsum)
  return(list(y = fit$time[fit$n.event==1], CIF = CIF))
}


### simulate competing risk data ###

simTwoCauseExpCSH <- function(n, v, h01, h02, beta1, beta2, 
                              c_scale = 1, mu = 0, X = NULL,
                              independent_c = TRUE,
                              censor = 0, seed = 1234){
  set.seed(seed)  # for reproducibility
  U <- runif(n) 
  if (is.null(X)) {
    X <- matrix(rnorm(n=n*length(beta1), mean = mu), nrow = n)
  }
  
  h <- h01*exp(X %*% beta1) + h02*exp(X %*% beta2)
  event.times <- (-(log(U)) / h)^(1/v)
  f.event<- apply(h01*exp(X %*% beta1)/h, 1, 
                  function(x) rbinom(n = 1, size = 1, prob = x))
  f.event<- ifelse(f.event== 0, 2, 1)
  
  if(independent_c==TRUE){
    cens.times <- rnorm(n, mean = 0, sd = sd(log(event.times)) * c_scale)
  }else{
    temp <- apply(X, 1, sum)
    cens.times <- sapply(temp, function(x){
      rnorm(n = 1, mean = x, sd = sd(log(event.times)) * c_scale)})
  }
  #cens.times <- log(runif(n, 0, 10) * c_scale)
  
  if(censor>0){
    result <- uniroot(find_mu_c, 
                      interval = c(-1e5, 1e5), 
                      censor = censor,
                      cens.times = cens.times, 
                      event.times = log(event.times), 
                      f.event = f.event, 
                      n = n)
    cens.times <- cens.times + result$root
    obs.times <- pmin(log(event.times), cens.times) %>% exp 
    obs.event <- c(log(event.times) <= cens.times) * f.event
  }else{
    obs.times <- event.times 
    obs.event <- f.event
  } 
  
  final_data<-data.frame(obs.times, obs.event, X) #%>% filter(obs.times>0)
  return(final_data)
}


find_mu_c <- function(c_mu, censor, cens.times, event.times, f.event, n) {
  cens.times <- cens.times + c_mu
  obs.times <- round(pmin(event.times, cens.times), 7)
  obs.event <- (event.times <= cens.times) * f.event
  result <- mean(obs.event == 0) - censor
  if (length(result) > 1) {
    warning("Result is not a scalar, taking mean.")
    result <- mean(result)
  }
  return(result)
}


simulateTwoCauseFineGrayModel <- function (n, v, beta1, beta2, lambda1 = 1, 
                                           X = NULL, mu = 0, p = 0.7,
                                           c_scale = 1, censor = 0, 
                                           independent_c = TRUE,
                                           seed = 1234) {
  set.seed(seed)
  if (length(beta1) != length(beta2)) 
    stop("Dimension of beta1 and beta2 should be the same")
  ncovs <- length(beta1)
  if (is.null(X)) {
    X <- matrix(rnorm(n * ncovs, mean = mu), nrow = n)
  }
  c.ind <- 1 + rbinom(n, 1, prob = (1 - p)^exp(X %*% beta1))
  ftime <- numeric(n)
  eta1 <- X[c.ind == 1, ] %*% beta1
  eta2 <- X[c.ind == 2, ] %*% beta2
  u1 <- runif(length(eta1))
  t1 <- (-log(1 - (1 - (1 - u1 * (1 - (1 - p)^exp(eta1)))^
                     (1/exp(eta1)))/p)/lambda1)^(1/v)
  t2 <- rweibull(length(eta2), shape = v, scale = exp(eta2)^(-1/v))
  
  if(independent_c==TRUE){
    cens.times <- rnorm(n, mean = 0, sd = sd(log(t1)) * c_scale)
  }else{
    temp <- apply(X, 1, sum)
    cens.times <- sapply(temp, function(x){
      rnorm(n = 1, mean = x, sd = sd(log(t1)) * c_scale)})
  }
  
  ftime[c.ind == 1] <- t1
  ftime[c.ind == 2] <- t2
  
  if(censor>0){
    result <- uniroot(find_mu_c, 
                      interval = c(-1e5, 1e5), 
                      censor = censor,
                      cens.times = cens.times, 
                      event.times = log(ftime), 
                      f.event = c.ind, 
                      n = n)
    cens.times <- cens.times + result$root
    obs.times <- pmin(log(ftime), cens.times) %>% exp 
    obs.event <- c(log(ftime) <= cens.times) * c.ind
  }else{
    obs.times <- ftime 
    obs.event <- c.ind
  } 
  
  final_data<-data.frame(obs.times, obs.event, X)
  return(final_data)
}


# Define two simulation function
sim_csh_check <- function(final_data, tau = NULL){
  if(is.null(tau)) tau <- max(final_data$obs.times)
  
  cox.model1 <- tryCatch({coxph(Surv(obs.times,obs.event==1) ~ ., 
                                data = final_data)}, warning = handle_warning) 
  cox.model2 <- tryCatch({coxph(Surv(obs.times,obs.event==2) ~ .,
                                data = final_data)}, warning = handle_warning) 
  if(is.null(cox.model1) | is.null(cox.model2)) return()
  
  model.info1 <- get_CIF(final_data$obs.times, final_data$obs.event, 
                         cox.model1, list(cox.model2), final_data)
  if(is.null(model.info1)) return()
  
  pred <- apply(model.info1$CIF, 2, m_cif, 
                time.cif = model.info1$y, 
                tau = tau)
  PAmeas <- pam.censor.cr(final_data$obs.times, final_data$obs.event,
                          tau, model.info1$CIF, model.info1$y)
  C_ind <- CindexCR(time = final_data$obs.times, 
                    status = final_data$obs.event,
                    predicted = pred, Cause_int = 1) %>% round(4)
  
  return(tibble(R.squared = PAmeas$R.squared, L.squared = PAmeas$L.squared,
                Psuedo.R = PAmeas$Psuedo.R, C_ind, model = "CSH"))
}

handle_warning <- function(w) {
  print(w)
  return()
}

handle_error <- function(e) {
  print(e)
  return()
}

sim_fng_check <- function(final_data, tau = NULL){
  if(is.null(tau)) tau <- max(final_data$obs.times)
  
  z <- tryCatch({crr(final_data$obs.times, 
                     final_data$obs.event, 
                     final_data[,-c(1,2)])}, error = handle_error)
  if(is.null(z)) return(NULL)
  if(z$converged==FALSE) return(NULL)
  
  z.p <- predict(z, as.matrix(final_data[,-c(1,2)]))
  if(is.null(dim(z.p[, -1]))) return(NULL)
  
  pred <- apply(z.p[, -1], 2, m_cif, time.cif = z.p[, 1], tau = tau)
  PAmeas <- pam.censor.cr(final_data$obs.times, final_data$obs.event,
                          tau, z.p[, -1], z.p[, 1])
  C_ind <- CindexCR(time = final_data$obs.times, 
                    status = final_data$obs.event,
                    predicted = pred, Cause_int = 1) %>% round(4)
  return(tibble(R.squared = PAmeas$R.squared, L.squared = PAmeas$L.squared,
                Psuedo.R = PAmeas$Psuedo.R, C_ind, model = "FnG"))
}

sim_function <- function(n, v, h01 = NULL, h02 = NULL, 
                         beta1, beta2, censor, p = NULL,
                         mu= -1, tau = NULL, tau_q = NULL,
                         independent_c = TRUE,
                         CSH = TRUE, seed = 1234){
  if(CSH == TRUE){
    final_data <- simTwoCauseExpCSH(n, v, h01, h02, beta1, beta2,
                                    independent_c = independent_c, 
                                    censor = censor, seed = seed, mu = mu)
  }else{
    final_data <- simulateTwoCauseFineGrayModel(n, v, beta1, beta2, p, mu = mu,
                                                independent_c = independent_c,
                                                censor = censor, seed = seed)
  }
  
  if(mean(final_data$obs.event==1)<0.001) return()
  out <- rbind(sim_csh_check(final_data, tau), sim_fng_check(final_data, tau))
  out$censor <- censor
  out$v <- v
  out$tau <- tau_q
  out$true_model <- ifelse(CSH == TRUE, "CSH", "FnG")
  out$beta1 <- paste0(beta1, collapse = ", ")
  out$beta2 <- paste0(beta2, collapse = ", ")
  out$h01 <- h01
  out$h02 <- h02
  out$n <- n
  out$p <- p
  return(out)
}

sim_quantile_t_fng <- function(param){
  final_data_fng <- simulateTwoCauseFineGrayModel(n = 5000, 
                                                  v = unlist(param[5]), 
                                                  beta1 = unlist(param[1]),
                                                  beta2 = unlist(param[2]),
                                                  p = unlist(param[4]), 
                                                  mu = 0, censor = 0, 
                                                  seed = unlist(param[6]))
  
  tibble(matrix(quantile(final_data_fng$obs.times,
                         c(0.25, 0.5, 0.75, 0.95, 1)), nrow = 1), 
         param, model = "fng", censor = mean(final_data_fng$obs.event==1))
}

sim_quantile_t_csh <- function(param){
  final_data_csh <- simTwoCauseExpCSH(n = 5000, 
                                      v = unlist(param[5]), 
                                      beta1 = unlist(param[1]),
                                      beta2 = unlist(param[2]),
                                      h01 = unlist(param[4])[1], 
                                      h02 = unlist(param[4])[2],
                                      mu = 0, censor = 0, 
                                      seed = unlist(param[6]))
  
  tibble(matrix(quantile(final_data_csh$obs.times,
                         c(0.25, 0.5, 0.75, 0.95, 1)), nrow = 1), 
         param, model = "csh", censor = mean(final_data_csh$obs.event==1))
}

generate_plot_df <- function(final_data, CSH=TRUE, tau = NULL){
  if(is.null(tau)) tau <- max(final_data$obs.times)
  if(CSH==TRUE){
    R.square <- sim_csh_check(final_data, tau = tau)
    
    cox.model1 <- coxph(Surv(obs.times,obs.event==1) ~ ., data = final_data) 
    
    cox.model2 <- coxph(Surv(obs.times,obs.event==2) ~ ., data = final_data)
    
    model.info1 <- get_CIF(final_data$obs.times, final_data$obs.event, 
                           cox.model1, list(cox.model2), final_data)
    
    pred1 <- apply(model.info1$CIF, 2, m_cif, 
                   time.cif = model.info1$y, tau = tau)
    rs1 <- as.matrix(final_data[,-c(1,2)]) %*% cox.model1$coefficients
    
    Cindex1 <- C_cr(time = final_data$obs.times, 
                    status = final_data$obs.event,
                    predicted = pred1, Time = tau,
                    Cause_int = 1) %>% round(4)
    
    df_plot <- data.frame(time = final_data$obs.times, 
                          cause = final_data$obs.event, 
                          pred1, rs1) %>% 
      pivot_longer(
        cols = c(pred1, rs1),
        names_to = c(".value", "group"),
        names_pattern = "(pred|rs)(\\d)"
      ) %>% 
      mutate(group =  paste0("R square: ", R.square$R.squared, 
                             ", C index: ", Cindex1),
             tau = tau)
    return(df_plot)
  } else {
    R.square <- sim_fng_check(final_data, tau = tau)
    z <- crr(final_data$obs.times, final_data$obs.event, final_data[,-c(1,2)])
    z.p <- predict(z, as.matrix(final_data[,-c(1,2)]))
    pred2 <- apply(z.p[,-1], 2, m_cif, 
                   time.cif = z.p[,1], tau = tau)
    rs2 <- as.matrix(final_data[,-c(1,2)]) %*% z$coef
    Cindex2 <- C_cr(time = final_data$obs.times, 
                    status = final_data$obs.event, 
                    predicted = pred2, Time = tau,
                    Cause_int = 1) %>% round(4)
    df_plot <- data.frame(time = final_data$obs.times, 
                          cause = final_data$obs.event, 
                          pred2, rs2) %>% 
      pivot_longer(
        cols = c(pred2, rs2),
        names_to = c(".value", "group"),
        names_pattern = "(pred|rs)(\\d)"
      ) %>% 
      mutate(group =  paste0("R square: ", R.square$R.squared, 
                             ", C index: ", Cindex2),
             tau = tau)
  }
  
}


#############################################################################
### C index competing risk sourced from SurvMetrics  (modified here) ########
#############################################################################
C_cr <- function (time, status, predicted, Time = NULL, Cause_int = 1) {
  if (any(is.na(time))) {
    stop("The input vector cannot have NA")
  }
  if (any(is.na(status))) {
    stop("The input vector cannot have NA")
  }
  if (any(!(status %in% c(0, 1, 2)))) {
    stop("The status must be 0 or 1 or 2")
  }
  if (any(is.na(predicted))) {
    stop("The input vector cannot have NA")
  }
  if (!(Cause_int %in% status)) {
    stop("Invalid input of Cause_int")
  }
  if (min(time) <= 0) {
    stop("Survival time must be positive")
  }
  Time_survival <- time
  Censoring <- ifelse(status == 0, 0, 1)
  Cause <- ifelse(status == 2, 2, 1)
  Prediction <- -log(predicted)
  if(is.null(Time)) Time <- max(Time_survival) + 1
  n <- length(Prediction)
  A <- matrix(0, nrow = n, ncol = n)
  B <- matrix(0, nrow = n, ncol = n)
  Q <- matrix(0, nrow = n, ncol = n)
  N_t <- matrix(0, nrow = n, ncol = n)
  Num_mat <- matrix(0, nrow = n, ncol = n)
  Den_mat <- matrix(0, nrow = n, ncol = n)
  Num <- 0
  Den <- 0
  for (i in 1:n) {
    A[i, which(Time_survival[i] < Time_survival)] <- 1
    B[i, intersect(intersect(which((Time_survival[i] >= 
                                      Time_survival)), which(Cause != Cause_int)), 
                   which(Censoring == 1))] <- 1
    Q[i, which(Prediction[i] > Prediction)] <- 1
  }
  for (i in 1:n) {
    if (Time_survival[i] <= Time && Cause[i] == Cause_int && 
        Censoring[i] == 1) {
      N_t[i, ] <- 1
    }
  }
  Num_mat <- (A + B) * Q * N_t
  Den_mat <- (A + B) * N_t
  Num <- sum(Num_mat)
  Den <- sum(Den_mat)
  return(Num/Den)
}