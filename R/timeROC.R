Compute.iid.KM <- function(times,status){
  #browser()
  times <- times[order(times)]
  status <- status[order(times)] 
  n <- length(times)
  mat.data<-cbind(times,as.numeric(status==0))
  colnames(mat.data)<-c("T","indic.Cens")
  # compute the empirical survival function corresponding to the counting process 1(\tilde{eta}=0, \tilde{T}<=t)
  hatSdeltaCensTc<-1-cumsum(mat.data[,c("indic.Cens")])/n  
  # Build the matrix required for computing  dM_C(u) for all time u (all observed times \tilde{T}_i)
  temp1 <- cbind(mat.data[,c("T","indic.Cens")],1-(1:n)/n,hatSdeltaCensTc)
  temp1 <- rbind(c(0,0,1,1),temp1) # Add the first row corresponding to time t=0
  colnames(temp1)<-c("T","indic.Cens","hatSTc","hatSdeltaCensTc")
  # compute hazard function of the censoring
  lambdaC<-(temp1[-1,"indic.Cens"])/(n:1)  
  # Add the column of the hazard function of the censoring (equal to 0 at time t=0)
  temp1<-cbind(temp1,c(0,lambdaC))
  colnames(temp1)[ncol(temp1)]<-"lambdaC"
  # Cumulative hazard of censoring
  LambdaC<-cumsum(lambdaC)         
  # Add the column of the cumulative hazard function of the censoring (equal to 0 at time t=0)
  temp1 <- cbind(temp1,c(0,LambdaC))
  colnames(temp1)[ncol(temp1)]<-"LambdaC"
  temp2<-temp1[-1,]
  # compute  martingale of censoring \hat{M}_{C_i}(u) for all time u (all observed times \tilde{T}_i) using previous matrix
  # We obtain a matrix. Each column contains the vector of M_{C_i}(\tilde{T}_j) for  all j.
  hatMC<-matrix(NA,n,n)
  for (i in 1:n){
    hatMC[,i] <-temp2[i,2]*as.numeric(temp2[i,1]<=temp2[,"T"])- c(temp2[0:i,"LambdaC"], rep(temp2[i,6],(n-i)))
  }  
  # In order to draw martingale paths
  #matplot(mat.data[,"T"],hatMC,type="l")
  #lines(mat.data[,"T"],rowMeans(hatMC),lwd=5)  
  # Compute d \hat{M}_{C_i} (u) for all time u (all observed times \tilde{T}_i)
  dhatMC<-rbind(hatMC[1,],hatMC[-1,]-hatMC[-nrow(hatMC),])
  # Compute d \hat{M}_{C_i} (u)/(S_{\tilde{T}}(u)) for all time u (all observed times \tilde{T}_i)
  # We need this for integrals in the martingale representation of the Kaplan-Meier estimator of the censoring survival function
  # function to divide d \hat{M}_{C_i} (u) by (S_{\tilde{T}}(u))
  MulhatSTc<-function(v){
    n <- length(v)
    v/c(1,1-(1:(n-1))/n)      # c(1,1-(1:(n-1))/n) is the at risk probability (S_{\tilde{T}}(u))
  }
  # apply the function for each column (corresponding to the
  # vector M_{C_i}(u)  for all time u (all observed times \tilde{T}_i), 
  # time \tilde{T}_i corresponds to the i-th row of the matrix)
  dhatMCdivST<-apply(dhatMC,2,MulhatSTc)
  # Compute \int_0^{\tilde{T}_j} d{ \hat{M}_{C_l} (u) } / (S_{\tilde{T}}(u)) for each subject l, we compute for all time \tilde{T}_j.
  # l=column, j=row
  MatInt0TcidhatMCksurEff<-apply(dhatMCdivST,2,cumsum)  # (Remark : on of the row corresponds to the previous step...) 
  colnames(MatInt0TcidhatMCksurEff)<-paste("M_{C_",1:length(times),"}",sep="")
  rownames(MatInt0TcidhatMCksurEff)<-times  
  return(MatInt0TcidhatMCksurEff)  
}

compute_iid_decomposition<-function(t,n,cause,F01t,St,weights,T,delta,marker,MatInt0TcidhatMCksurEff){ 
  # indicator vectors 
  Cases<-(T< t & delta==cause)
  Controls_1<-(T> t)
  Controls_2<-(T< t &  delta!=cause & delta!=0)
  if(sum(Controls_2)>0){
    compute_iid_decomposition_competing_risks(t,n,cause,F01t,St,weights,T,delta,marker,MatInt0TcidhatMCksurEff) 
  }else{
    compute_iid_decomposition_survival(t,n,cause,F01t,St,weights,T,delta,marker,MatInt0TcidhatMCksurEff)
  }  
}

compute_iid_decomposition_survival<-function(t,n,cause,F01t,St,weights,T,delta,marker,MatInt0TcidhatMCksurEff){  
  start_total<-Sys.time()
  # indicator vectors 
  Cases<-(T< t & delta==cause)
  Controls_1<-(T> t )
  # vectors which indicates the indexes of Cases and the Controls
  which_Cases<-which(T< t & delta==cause)
  which_Controls_1<-which(T> t )
  # compute the weights. 
  Weights_cases_all<-1/(weights$IPCW.subjectTimes*n)
  Weights_cases<-Weights_cases_all
  Weights_cases[!Cases]<-0  #(0 if not a case)
  Weights_controls_1<-rep(1/(weights$IPCW.times[which(weights$times==t)]*n),times=n)
  Weights_controls_1[!Controls_1]<-0  #(0 if not a control)
  # compute vector indicator of censoring (event is censoring !)
  indic_Cens<-as.numeric(delta==0)
  # compute the matrix with all information. The matrix is order by order(t)
  Mat_data<-cbind(T,delta,indic_Cens,marker,Cases,Controls_1,Weights_cases,Weights_controls_1)
  
  ## MatInt0TcidhatMCksurEff <- Compute.iid.KM(times=T,status=delta)
  # {{{ STEP : Compute terms  {\hat{h}_{tij}}_1 and {\hat{h}_{tij}}_2
  #start_htij<-Sys.time() 
  # function that eats the matrix W1 (defined just after) that depends on subject i and returns 
  # the vector of {\hat{h}_{tij}}_1 
  htij1<-function(V,tps=t){
    as.numeric(V[,1]>tps)*(as.numeric(V[,4]>V[,2]) + 0.5*as.numeric(V[,4]==V[,2])) *(V[,3]*V[,5])*(n*n)
  }
  # compute frequencies of cases and controls to define 
  #the size of the matrix  Mathtij1 
  nb_Cases<-sum(T< t & delta==cause)
  nb_Controls_1<-sum(T> t )
  # To save computation time, we loop only on control 1 for Mathtij1 
  Mat_data_cont1<-Mat_data[which_Controls_1,]
  # initialise  Mathtij1  with its right size !
  Mathtij1<-matrix(NA,nb_Controls_1,nb_Cases)
  # loop on all cases i. We loop only on Cases to save computation time !  
  for (i in which_Cases){
    W1<-cbind(Mat_data_cont1[,c("T","marker")],
              rep(Mat_data[i,c("Weights_cases")],nb_Controls_1),
              rep(Mat_data[i,c("marker")],nb_Controls_1),
              Mat_data_cont1[,c("Weights_controls_1")])
    # fill the column i of  Mathtij1 and  Mathtij2
    Mathtij1[,which(i==which_Cases)]<-htij1(W1) 
  }
  # matrix Mathtij1  : i for columns, j for rows
  #browser() # nice function for debugging !
  #stop_htij<-Sys.time()
  #print(difftime(stop_htij,start_htij,units="sec"))
  # compute \hat{h}_t
  ht<-(sum(Mathtij1) )/(n*n) 
  # vector of \hat{f}_{i1t}
  vect_dit<-as.numeric(Mat_data[,c("T")]<=t)*as.numeric(Mat_data[,c("delta")]==cause)*Mat_data[,c("Weights_cases")]*n
  # We can check we have F01t by mean(vect_dit)
  #print("F01t ??")
  #print(c(mean(vect_dit),F01t))
  # }}}
  # {{{ Final step : to compute iid representation of AUC^*(t)
  start_iid_AUC1<-Sys.time()
  # Compute the vecor of all sum_{i=1}^n of {\hat{h}_{tij}}_1 for all j
  colSums_Mathtij1<-rep(0,n) # initialise at 0
  colSums_Mathtij1[which_Cases]<-colSums(Mathtij1) # when i is a case,  then we sum the column of  Mathtij1  
  # Compute the vecor of all sum_{j=1}^n of {\hat{h}_{tij}}_1 for all i  
  rowSums_Mathtij1<-rep(0,n) # initialize at 0
  rowSums_Mathtij1[which_Controls_1]<-rowSums(Mathtij1)# when  j is a control 1, then we sum the row of  Mathtij1
  hathtstar<-(sum(Mathtij1)  )/(n*n)  
  #print("AUC1 ???")
  #print(hathtstar/(F01t*St))  
  # compute the vector of \frac{1_{\tilde{T}_i>=t}}{ \hat{S}_{\tilde{T}}(t)}
  vect_Tisupt<-as.numeric(Mat_data[,c("T")]>t)/( sum(as.numeric(Mat_data[,c("T")]>t))/n )   
  sum_ij_a_k_fixe<-function(k){
    Pour_sum_ij_a_k_fixe<- t(Mathtij1)*(1+MatInt0TcidhatMCksurEff[which_Cases,k]) 
    Pour_sum_ij_a_k_fixe_3<-vect_dit*(1+MatInt0TcidhatMCksurEff[,k])
    Pour_sum_ij_a_k_fixe_3b<-(hathtstar)*(  vect_Tisupt    +  (1/F01t)*(Pour_sum_ij_a_k_fixe_3-F01t) )
    La_sum_ij_a_k_fixe<- sum(Pour_sum_ij_a_k_fixe)/n - sum(Pour_sum_ij_a_k_fixe_3b) 
    return(La_sum_ij_a_k_fixe)
  }  
  #print("F01t*St")
  #print(F01t*St)  
  Les_sum_ij_a_k_fixe<-(sapply(1:n,sum_ij_a_k_fixe))/(F01t*St)  
  Les_sum_ik_a_j_fixe<-(rowSums_Mathtij1 - n*hathtstar)/(F01t*St)
  Les_sum_jk_a_i_fixe<- (colSums_Mathtij1 - n*hathtstar*(vect_Tisupt+(1/F01t)*(vect_dit-F01t)))/(F01t*St)
  # We compute the iid representation of the AUC estimator
  hatIFstar<- (Les_sum_ij_a_k_fixe + Les_sum_ik_a_j_fixe +  Les_sum_jk_a_i_fixe)/(n)
  stop_iid_AUC1<-Sys.time()
  # }}}
  # we compute the standard error of the AUC estimator
  seAUCstar<-sd(hatIFstar)/sqrt(n)
  #browser() # nice function for debugging
  stop_total<-Sys.time()
  total_time<-difftime(stop_total,start_total,units="secs")
  total_time_iid_AUC1<-difftime(stop_iid_AUC1,start_iid_AUC1,units="secs")
  computation_times<-c(total_time)
  names(computation_times)<-c("total_time")
  return(list(iid_representation_AUC=rep(NA,n),
              iid_representation_AUCstar=hatIFstar,
              seAUC=NA,seAUCstar=seAUCstar,
              computation_times=computation_times)
  )
}

compute_iid_decomposition_competing_risks<-function(t,n,cause,F01t,St,weights,T,delta,marker,MatInt0TcidhatMCksurEff){
  start_total<-Sys.time() 
  # indicator vectors 
  Cases<-(T< t & delta==cause)
  Controls_1<-(T> t )
  Controls_2<-(T< t &  delta!=cause & delta!=0)
  # vectors which indicates the indexes of Cases and the Controls
  which_Cases<-which(T< t & delta==cause)
  which_Controls_1<-which(T> t )
  which_Controls_2<-which(T< t &  delta!=cause & delta!=0)
  # compute the weights. 
  Weights_cases_all<-1/(weights$IPCW.subjectTimes*n)
  Weights_cases<-Weights_cases_all
  Weights_controls_2<-Weights_cases_all
  Weights_cases[!Cases]<-0  #(0 if not a case)
  Weights_controls_2[!Controls_2]<-0  #(0 if not a control)
  Weights_controls_1<-rep(1/(weights$IPCW.times[which(weights$times==t)]*n),times=n)
  Weights_controls_1[!Controls_1]<-0  #(0 if not a control) 
  # compute vector indicator of censoring (event is censoring !)
  indic_Cens<-as.numeric(delta==0)
  # compute the matrix with all information. The matrix is order by order(t)
  Mat_data<-cbind(T,delta,indic_Cens,marker,Cases,Controls_1,Controls_2,Weights_cases,Weights_controls_2,Weights_controls_1)
  
  ## MatInt0TcidhatMCksurEff <- Compute.iid.KM(times=T,status=delta)
  Int0tdMCsurEffARisk <- MatInt0TcidhatMCksurEff[max(which(Mat_data[,"T"]<=t)),]   
  # {{{ Step : Compute terms  {\hat{h}_{tij}}_1 and {\hat{h}_{tij}}_2
  #start_htij<-Sys.time() 
  # function that eats the matrix W1 (defined just after) that depends on subject i and returns 
  # the vector of {\hat{h}_{tij}}_1 
  htij1<-function(V,tps=t){
    as.numeric(V[,1]>tps)*(as.numeric(V[,4]>V[,2]) + 0.5*as.numeric(V[,4]==V[,2])) *(V[,3]*V[,5])*(n*n)
  }
  # function that eats the matrix W2 (defined just after) that depends on subject i and returns 
  # the vector of {\hat{h}_{tij}}_2
  htij2<-function(V,tps=t){
    as.numeric(V[,1]<=tps)*(as.numeric(V[,4]>V[,2]) + 0.5*as.numeric(V[,4]==V[,2]))*as.numeric(V[,6]!=0)*as.numeric(V[,6]!=cause) *(V[,3]*V[,5])*(n*n)
  }
  # compute frequencies of cases and controls to define 
  #the size of the matrix  Mathtij1 and  Mathtij1
  nb_Cases<-sum(T< t & delta==cause)
  nb_Controls_1<-sum(T> t )
  nb_Controls_2<-sum(T< t &  delta!=cause & delta!=0)  
  # To save computation time, we loop only on control 1 for Mathtij1 and 
  # only on control 2 for Mathtij2
  Mat_data_cont1<-Mat_data[which_Controls_1,]
  Mat_data_cont2<-Mat_data[which_Controls_2,]
  # initialise  Mathtij1 and  Mathtij2 with their right sizes !
  Mathtij1<-matrix(NA,nb_Controls_1,nb_Cases)
  Mathtij2<-matrix(NA,nb_Controls_2,nb_Cases) 
  # loop on all cases i. We loop only on Cases to save computation time !
  for (i in which_Cases){
    W1<-cbind(Mat_data_cont1[,c("T","marker")],
              rep(Mat_data[i,c("Weights_cases")],nb_Controls_1),
              rep(Mat_data[i,c("marker")],nb_Controls_1),
              Mat_data_cont1[,c("Weights_controls_1")])
    W2<-cbind(Mat_data_cont2[,c("T","marker")],
              rep(Mat_data[i,c("Weights_cases")],nb_Controls_2),
              rep(Mat_data[i,c("marker")],nb_Controls_2),
              Mat_data_cont2[,c("Weights_controls_2")],Mat_data_cont2[,c("delta")])
    # fill the column i of  Mathtij1 and  Mathtij2
    Mathtij1[,which(i==which_Cases)]<-htij1(W1) 
    Mathtij2[,which(i==which_Cases)]<-htij2(W2)
  }
  # matrix Mathtij1 and  Mathtij2 : i for columns, j for rows
  #browser() # nice function for debugging !
  #stop_htij<-Sys.time()
  #print(difftime(stop_htij,start_htij,units="sec"))
  # compute \hat{h}_t
  ht<-(sum(Mathtij1) +sum(Mathtij2) )/(n*n) 
  #print("ht") 
  #print(ht)
  # We can check we have the AUC by \hat{h}_t/((1-F01t)*F01t)
  #AUChtij<-ht/((1-F01t)*F01t)
  #print("check_AUC")
  #print(AUChtij)
  # vector of \hat{f}_{i1t}
  vect_dit<-as.numeric(Mat_data[,c("T")]<=t)*as.numeric(Mat_data[,c("delta")]==cause)*Mat_data[,c("Weights_cases")]*n 
  # We can check we have F01t by mean(vect_dit)
  #print("F01t ??")
  #print(c(mean(vect_dit),F01t))
  # }}} 
  # {{{ FINAL step : compute iid representation of AUC(t)
  # we compute this step only in presence of competing risks 
  start_iid_AUC2<-Sys.time()
  # Let' recall :
  # Mathtij1 # matrix of  {\hat{h}_{tij}}_1, i for columns, j for rows
  # Mathtij2 # matrix of  {\hat{h}_{tij}}_2, i for columns, j for rows
  # MatInt0TcidhatMCksurEff # a matrix of \int_0^{\tilde{T}_j} d{ \hat{M}_{C_l} (u) } / (S_{\tilde{T}}(u)),  l=column, j=row for
  # A function that eats index l and 
  # returns \frac{1}{n}\sum_{i=1}^n \sum_{j=1}^n \sum_{k=1}^n \Psi_{ijkl}(t)
  sum_ijk_a_l_fixe<-function(l){
    Pr_sum_ijk_a_l_fixe_1<-Mathtij1*(1+Int0tdMCsurEffARisk[l])
    Pr_sum_ijk_a_l_fixe_2<-Mathtij2* (1+MatInt0TcidhatMCksurEff[which_Controls_2,l])
    La_sum_ijk_a_l_fixe<- (sum(Pr_sum_ijk_a_l_fixe_1) + sum(Pr_sum_ijk_a_l_fixe_2)- n^2*ht)
    return(La_sum_ijk_a_l_fixe)
  }
  # A function that eats index k and 
  #returns \frac{1}{n}\sum_{i=1}^n \sum_{j=1}^n \sum_{l=1}^n \Psi_{ijkl}(t)
  sum_ijl_a_k_fixe<-function(k){
    Pour_sum_ijl_a_k_fixe_1<- t(Mathtij1)*(1+MatInt0TcidhatMCksurEff[which_Cases,k]) 
    Pour_sum_ijl_a_k_fixe_2<- t(Mathtij2)*(1+MatInt0TcidhatMCksurEff[which_Cases,k])   
    Pour_sum_ijl_a_k_fixe_3<-vect_dit*(1+MatInt0TcidhatMCksurEff[,k])
    Pour_sum_ijl_a_k_fixe_3b<-(ht*(1-2*F01t)/(F01t*(1-F01t)))*(Pour_sum_ijl_a_k_fixe_3-F01t)
    La_sum_ijl_a_k_fixe<-( (sum(Pour_sum_ijl_a_k_fixe_1) +sum(Pour_sum_ijl_a_k_fixe_2) )- n^2*ht -n*sum(Pour_sum_ijl_a_k_fixe_3b) )
    return(La_sum_ijl_a_k_fixe)
  }
  # Compute the vecor of all sum_{i=1}^n of {\hat{h}_{tij}}_1 for all j
  colSums_Mathtij1<-rep(0,n) # initialise at 0
  colSums_Mathtij1[which_Cases]<-colSums(Mathtij1) # when i is a case,  then we sum the column of  Mathtij1
  # Compute the vecor of all sum_{i=1}^n of {\hat{h}_{tij}}_2 for all j
  colSums_Mathtij2<-rep(0,n) # initialise at 0
  colSums_Mathtij2[which_Cases]<-colSums(Mathtij2) # when i is a case, then we sum the column of  Mathtij2
  # Compute the vecor of all sum_{j=1}^n of {\hat{h}_{tij}}_1 for all i  
  rowSums_Mathtij1<-rep(0,n) # initialize at 0
  rowSums_Mathtij1[which_Controls_1]<-rowSums(Mathtij1)# when  j is a control 1, then we sum the row of  Mathtij1
  # Compute the vecor of all sum_{j=1}^n of {\hat{h}_{tij}}_2 for all i
  rowSums_Mathtij2<-rep(0,n) # initialize at 0
  rowSums_Mathtij2[which_Controls_2]<-rowSums(Mathtij2) # when  j is a control 2, then we sum the row of  Mathtij2
  # we compute \frac{1}{n}\sum_{j=1}^n \sum_{k=1}^n \sum_{l=1}^n \Psi_{ijkl}(t)
  Les_sum_jkl_a_i_fixe<-( (colSums_Mathtij1 + colSums_Mathtij2)*n - n^2*ht - ( ht*n^2*(1-2*F01t) / (F01t*(1-F01t)) ) *(vect_dit - F01t) )/(F01t*(1-F01t))
  # we compute \frac{1}{n}\sum_{i=1}^n \sum_{k=1}^n \sum_{l=1}^n \Psi_{ijkl}(t)
  Les_sum_ikl_a_j_fixe<-((rowSums_Mathtij1 + rowSums_Mathtij2)*n - n^2*ht)/(F01t*(1-F01t))
  # we compute \frac{1}{n}\sum_{i=1}^n \sum_{j=1}^n \sum_{k=1}^n \Psi_{ijkl}(t) 
  Les_sum_ijk_a_l_fixe<-(sapply(1:n,sum_ijk_a_l_fixe))/(F01t*(1-F01t))
  #start_step<-Sys.time()
  # we compute \frac{1}{n}\sum_{i=1}^n \sum_{j=1}^n \sum_{l=1}^n \Psi_{ijkl}(t)
  Les_sum_ijl_a_k_fixe<-(sapply(1:n,sum_ijl_a_k_fixe))/(F01t*(1-F01t))
  #stop_step<-Sys.time()
  #print(difftime(stop_step,start_step,units="sec")) 
  # We compute the iid representation of the AUC estimator
  hatIF<- (Les_sum_jkl_a_i_fixe + Les_sum_ikl_a_j_fixe + Les_sum_ijk_a_l_fixe + Les_sum_ijl_a_k_fixe)/(n*n)
  stop_iid_AUC2<-Sys.time()
  # }}}
  # {{{ Step : compute iid representation of AUC^*(t) 
  start_iid_AUC1<-Sys.time()  
  hathtstar<-(sum(Mathtij1)  )/(n*n)
  #print("AUC1 ???")
  #print(hathtstar/(F01t*St))
  # compute the vector of \frac{1_{\tilde{T}_i>=t}}{ \hat{S}_{\tilde{T}}(t)}
  vect_Tisupt<-as.numeric(Mat_data[,c("T")]>t)/( sum(as.numeric(Mat_data[,c("T")]>t))/n ) 
  sum_ij_a_k_fixe<-function(k){
    Pour_sum_ij_a_k_fixe<- t(Mathtij1)*(1+MatInt0TcidhatMCksurEff[which_Cases,k]) 
    Pour_sum_ij_a_k_fixe_3<-vect_dit*(1+MatInt0TcidhatMCksurEff[,k])
    Pour_sum_ij_a_k_fixe_3b<-(hathtstar)*(  vect_Tisupt    +  (1/F01t)*(Pour_sum_ij_a_k_fixe_3-F01t) )
    La_sum_ij_a_k_fixe<- sum(Pour_sum_ij_a_k_fixe)/n - sum(Pour_sum_ij_a_k_fixe_3b) 
    return(La_sum_ij_a_k_fixe)
  }
  #print("F01t*St")
  #print(F01t*St)
  Les_sum_ij_a_k_fixe<-(sapply(1:n,sum_ij_a_k_fixe))/(F01t*St)
  Les_sum_ik_a_j_fixe<-(rowSums_Mathtij1 - n*hathtstar)/(F01t*St)
  Les_sum_jk_a_i_fixe<- (colSums_Mathtij1 - n*hathtstar*(vect_Tisupt+(1/F01t)*(vect_dit-F01t)))/(F01t*St)
  # We compute the iid representation of the AUC estimator
  hatIFstar<- (Les_sum_ij_a_k_fixe + Les_sum_ik_a_j_fixe +  Les_sum_jk_a_i_fixe)/(n)
  stop_iid_AUC1<-Sys.time()
  # }}}
  # we compute the standard error of the AUC estimators
  seAUC<-sd(hatIF)/sqrt(n)
  seAUCstar<-sd(hatIFstar)/sqrt(n)
  #browser() # nice function for debugging
  stop_total<-Sys.time()
  total_time<-difftime(stop_total,start_total,units="secs")
  total_time_iid_AUC1<-difftime(stop_iid_AUC1,start_iid_AUC1,units="secs")
  total_time_iid_AUC2<-difftime(stop_iid_AUC2,start_iid_AUC2,units="secs")
  additional_times<-c(total_time_iid_AUC1,total_time_iid_AUC2)
  computation_times<-c(total_time)
  names(computation_times)<-c("total_time")
  return(list(iid_representation_AUC=hatIF,
              iid_representation_AUCstar=hatIFstar,
              seAUC=seAUC,seAUCstar=seAUCstar,
              computation_times=computation_times)
  )
}

timeROC<-function(T,delta,marker,other_markers=NULL,cause,weighting="marginal",times,ROC=TRUE,iid=FALSE){
  # {{{ check some inputs
  if (length(delta)!=length(T) | length(marker)!=length(T) | length(delta)!=length(T)){
    stop("lengths of vector T, delta and marker have to be equal\n") }
  if (missing(times)){
    stop("Choose at least one time for computing the time-dependent AUC\n") } 
  if (!weighting %in% c("marginal","cox","aalen")){
    stop("the weighting argument must be marginal (default), cox or aalen.\n") }  
  if (weighting %in% c("cox","aalen") & !missing(other_markers) & !("matrix" %in% class(other_markers))){
    stop("argument other_markers must be a matrix\n") }
  if (weighting %in% c("cox","aalen") & !missing(other_markers)){
    if(!nrow(other_markers)==length(marker))  stop("lengths of vector T, delta, marker and number of rows of other_markers have to be equal\n")
  }
  # }}}
  # {{{ check if there are missing values, and delete rows with missing values
  if (weighting %in% c("cox","aalen") & !missing(other_markers) ){
    is_not_na<-as.logical(apply(!is.na(cbind(T,delta,marker,other_markers)),1,prod))
    T<-T[is_not_na]
    delta<-delta[is_not_na]
    marker<-marker[is_not_na]
    other_markers<-as.matrix(other_markers[is_not_na,])
  }else{
    is_not_na<-as.logical(apply(!is.na(cbind(T,delta,marker)),1,prod)) 
    T<-T[is_not_na]
    delta<-delta[is_not_na]
    marker<-marker[is_not_na]
  }
  # }}} 
  start_computation_time<-Sys.time()
  # {{{ create some usefull objects
  n<-length(T)
  n_marker<-length(unique(marker))
  n_times<-length(times)
  if (n_times==1){times<-c(0,times)
  n_times<-2}           # trick to use ipcw.cox() even if there is only one time
  times<-times[order(times)]
  times_names<-paste("t=",times,sep="")
  # }}}
  # {{{ output initialisation
  AUC_1<-rep(NA,n_times)
  AUC_2<-rep(NA,n_times)
  CumInci<-rep(NA,n_times)
  surv<-rep(NA,n_times)
  names(AUC_1)<-times_names
  names(AUC_2)<-times_names
  names(CumInci)<-times_names
  names(surv)<-times_names
  Stats<-matrix(NA,nrow=n_times,ncol=4)
  colnames(Stats)<-c("Cases","survivor at t","Other events at t","Censored at t")
  rownames(Stats)<-times_names
  # }}}
  # {{{  computation of weights (1/2)
  # we need to order to use the pec::ipcw() fonction
  order_T<-order(T)
  T <- T[order_T]
  delta <- delta[order_T]
  marker<- marker[order_T]
  # use ipcw function from pec package
  if(weighting=="marginal"){
    weights <- pec::ipcw(Surv(failure_time,status)~1,data=data.frame(failure_time=T,status=as.numeric(delta!=0)),method="marginal",times=times,subjectTimes=T,subjectTimesLag=1)
  }
  if(weighting=="cox"){
    if (missing(other_markers)){marker_censoring<-marker } 
    other_markers<-other_markers[order_T,]
    marker_censoring<-cbind(marker,other_markers)
    colnames(marker_censoring)<-paste("X", 1:ncol(marker_censoring), sep="")
    fmla <- as.formula(paste("Surv(T,status)~", paste(paste("X", 1:ncol(marker_censoring), sep=""), collapse= "+")))
    data_weight<-as.data.frame(cbind(data.frame(T=T,status=as.numeric(delta!=0)),marker_censoring))
    weights <- pec::ipcw(fmla,data=data_weight,method="cox",times=as.matrix(times),subjectTimes=data_weight[,"T"],subjectTimesLag=1)
  }
  if(weighting=="aalen"){
    if (missing(other_markers)){marker_censoring<-marker }
    other_markers<-other_markers[order_T,]
    marker_censoring<-cbind(marker,other_markers)
    colnames(marker_censoring)<-paste("X", 1:ncol(marker_censoring), sep="")
    fmla <- as.formula(paste("Surv(T,status)~", paste(paste("X", 1:ncol(marker_censoring), sep=""), collapse= "+")))
    data_weight<-as.data.frame(cbind(data.frame(T=T,status=as.numeric(delta!=0)),marker_censoring))
    weights <- pec::ipcw(fmla,data=data_weight,method="aalen",times=as.matrix(times),subjectTimes=data_weight[,"T"],subjectTimesLag=1)
  }
  # we order by marker values (in order to compute Se and Sp)
  order_marker<-order(-marker)
  Mat_data<-cbind(T,delta,marker)[order_marker,]
  colnames(Mat_data)<-c("T","delta","marker")
  # Create some weights
  Weights_cases_all<-1/(weights$IPCW.subjectTimes*n)
  Weights_cases_all<-Weights_cases_all[order_marker]
  # }}}
  # {{{ Make TP and FP outputs if needed
  if(ROC==TRUE){ 
    FP_1<-matrix(NA,nrow=(n_marker+1),ncol=n_times)
    TP<-matrix(NA,nrow=(n_marker+1),ncol=n_times)
    FP_2<-matrix(NA,nrow=(n_marker+1),ncol=n_times)
    colnames(FP_1)<-times_names
    colnames(TP)<-times_names
    colnames(FP_2)<-times_names
  } else{FP_1<-NA
  FP_2<-NA
  TP<-NA}
  # }}}
  # {{{ loop on all timepoints t
  for(t in 1:n_times){
    Cases<-(Mat_data[,"T"]< times[t] &  Mat_data[,"delta"]==cause)
    Controls_1<-(Mat_data[,"T"]> times[t] )
    Controls_2<-(Mat_data[,"T"]< times[t] &  Mat_data[,"delta"]!=cause & Mat_data[,"delta"]!=0)  
    if (weights$method!="marginal"){ 
      Weights_controls_1<-1/(weights$IPCW.times[,t]*n)  }
    else{
      Weights_controls_1<-rep(1/(weights$IPCW.times[t]*n),times=n)
    }
    Weights_controls_1<-Weights_controls_1[order_marker] 
    Weights_cases<-Weights_cases_all 
    Weights_controls_2<-Weights_cases_all
    Weights_cases[!Cases]<-0
    Weights_controls_1[!Controls_1]<-0
    Weights_controls_2[!Controls_2]<-0
    den_TP_t<-sum(Weights_cases)
    den_FP_1_t<-sum(Weights_controls_1)
    den_FP_2_t<-sum(Weights_controls_2)+sum(Weights_controls_1)
    if(den_TP_t!=0){  
      TP_tbis<-c(0,cumsum(Weights_cases))/den_TP_t
      TP_t<-TP_tbis[!duplicated(marker[order_marker])]
    }
    else TP_t<-NA
    if(den_FP_1_t!=0){
      FP_1_tbis<-c(0,cumsum(Weights_controls_1))/den_FP_1_t
      FP_1_t<-FP_1_tbis[!duplicated(marker[order_marker])]}
    else FP_1_t<-NA
    if(den_FP_2_t!=0){
      FP_2_tbis<-c(0,cumsum(Weights_controls_1)+cumsum(Weights_controls_2))/den_FP_2_t
      FP_2_t<-FP_2_tbis[!duplicated(marker[order_marker])]}
    else FP_2_t<-NA
    # internal fonction to compute an area under a curve by trapezoidal rule
    AireTrap<-function(Abs,Ord){
      nobs<-length(Abs)
      dAbs<-Abs[-1]-Abs[-nobs]
      mil<-(Ord[-nobs]+Ord[-1])/2
      area<-sum(dAbs*mil)
      return(area)
    }
    if ( den_TP_t*den_FP_1_t != 0){AUC_1[t]<-AireTrap(FP_1_t,TP_t)}
    else AUC_1[t]<-NA
    if ( den_TP_t*den_FP_2_t != 0){AUC_2[t]<-AireTrap(FP_2_t,TP_t)}
    else AUC_2[t]<-NA
    if(ROC==TRUE){ 
      TP[,t]<-TP_t
      FP_1[,t]<-FP_1_t
      FP_2[,t]<-FP_2_t
    }  
    CumInci[t]<-c(den_TP_t)
    surv[t]<-c(den_FP_1_t)
    Stats[t,]<-c(sum(Cases),sum(Controls_1),sum(Controls_2),n-sum(Cases)-sum(Controls_1)-sum(Controls_2))
  }
  # }}}
  inference<-NA
  if (iid==TRUE){   
    if(weighting!="marginal"){
      stop("Error : Weighting must be marginal for computing the iid representation \n Choose iid=FALSE or weighting=marginal in the input arguments")
    } 
    else{      
      # create iid representation required for inference procedures
      out_iid<-vector("list", n_times)
      names(out_iid)<-paste("t=",times,sep="")
      vect_iid_comp_time<-rep(NA,times=n_times)
      names(vect_iid_comp_time)<-paste("t=",times,sep="")
      mat_iid_rep<-matrix(NA,nrow=n,ncol=n_times)
      colnames(mat_iid_rep)<-paste("t=",times,sep="")  
      mat_iid_rep_star<-matrix(NA,nrow=n,ncol=n_times)
      colnames(mat_iid_rep_star)<-paste("t=",times,sep="")  
      vetc_se<-rep(NA,times=n_times)
      names(vetc_se)<-paste("t=",times,sep="")
      vetc_sestar<-rep(NA,times=n_times)    
      names(vetc_sestar)<-paste("t=",times,sep="")   
      # compute iid for Kaplan Meier
      MatInt0TcidhatMCksurEff <- Compute.iid.KM(times=T,status=delta)
      for (j in 1:n_times){
        #compute iid representation when AUC can be computed
        if(!is.na(AUC_1[j]) | !is.na(AUC_2[j])){
          out_iid[[j]]<-compute_iid_decomposition(t=times[j],n=n,cause=cause,F01t=CumInci[j],St=surv[j],weights,T,delta,marker,MatInt0TcidhatMCksurEff=MatInt0TcidhatMCksurEff)} 
        else{
          out_iid[[j]]<-NA}
        #browser()
        #save output for inference for AUC_1 when AUC_1 can be computed        
        if(!is.na(AUC_1[j])){
          mat_iid_rep_star[,j]<-out_iid[[j]]$iid_representation_AUCstar
          vetc_sestar[j]<-out_iid[[j]]$seAUCstar
          vect_iid_comp_time[j]<-out_iid[[j]]$computation_times               
        }
        #save output for inference for AUC_2 when AUC_2 can be computed         
        if(!is.na(AUC_2[j])){
          mat_iid_rep[,j]<-out_iid[[j]]$iid_representation_AUC
          vetc_se[j]<-out_iid[[j]]$seAUC
          vect_iid_comp_time[j]<-out_iid[[j]]$computation_times               
        }   
      }
      inference<-list(mat_iid_rep_2=mat_iid_rep,
                      mat_iid_rep_1=mat_iid_rep_star,
                      vect_sd_1=vetc_sestar,
                      vect_sd_2=vetc_se,
                      vect_iid_comp_time=vect_iid_comp_time
      )
    }
  }
  stop_computation_time<-Sys.time() 
  # output if there is competing risks or not
  if (max(Stats[,3])==0){
    out <- list(TP=TP,FP=FP_1,AUC=AUC_1,times=times,
                CumulativeIncidence=CumInci,survProb=surv,n=n,Stats=Stats[,c(1,2,4)],weights=weights,
                inference=inference,computation_time=difftime(stop_computation_time,start_computation_time,units="secs"),iid=iid)
    class(out) <- "ipcwsurvivalROC"
    out
  }else{
    out <- list(TP=TP,FP_1=FP_1,AUC_1=AUC_1,FP_2=FP_2,AUC_2=AUC_2,times=times,
                CumulativeIncidence=CumInci,survProb=surv,n=n,Stats=Stats,weights=weights,
                inference=inference,computation_time=difftime(stop_computation_time,start_computation_time,units="secs"),iid=iid)
    class(out) <- "ipcwcompetingrisksROC"
    out    
  }
}