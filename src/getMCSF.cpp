#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List getMCSF(Eigen::VectorXd & beta, Eigen::VectorXd & gamma1, 
                 Eigen::VectorXd & alpha1,
                 Eigen::MatrixXd & H01,
                 Eigen::MatrixXd & Sig, double sigma, 
                 const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
                 const Eigen::VectorXd & Y, const Eigen::MatrixXd & X2, 
                 const Eigen::VectorXd & survtime, 
                 const Eigen::VectorXd & cmprsk, const Eigen::VectorXd & mdata, 
                 const Eigen::VectorXi & mdataS,
                 const Eigen::MatrixXd & FUNBS, 
                 const Eigen::MatrixXd & FUNEC, 
                 const Eigen::MatrixXd & FUNBEC,
                 const Eigen::MatrixXd & FUNBSEC, 
                 const Eigen::MatrixXd & FUNB){
  
  int a =H01.rows();
  int k = mdataS.size();
  int p1 = X1.cols();
  int p1a = Z.cols();
  
  int i,q,j,t,u;
  
  double scalef=0;
  
  Eigen::MatrixXd X1X1T  = Eigen::MatrixXd::Zero(p1, p1);
  Eigen::VectorXd X11 = Eigen::VectorXd::Zero(p1);
  Eigen::MatrixXd bbT  = Eigen::MatrixXd::Zero(p1a, p1a);
  Eigen::MatrixXd bbT2  = Eigen::MatrixXd::Zero(p1a, p1a);
  
  /* calculate beta */
  for (j=0;j<k;j++) {
    
    q=mdata(j);
    
    for (i=0;i<q;i++) {
      
      X1X1T += MultVVoutprod(X1.row(mdataS(j)-1+i));
      X11 += (Y(mdataS(j)-1+i) - MultVV(Z.row(mdataS(j)-1+i), FUNB.col(j)))*X1.row(mdataS(j)-1+i);
                                                          
    }
    
  }
  
  beta = X1X1T.inverse()*X11;
  
  /* calculate Sig*/
  Sig = Eigen::MatrixXd::Zero(p1a, p1a);
  
  for (j=0;j<k;j++) {
    
    for(t=0;t<p1a;t++) Sig(t,t)+=FUNBS(t,j);
    if (p1a>1) {
      for(q=1;q<p1a;q++)
      {
        for(t=0;t<(p1a-q);t++) {
          Sig(t,q+t)+=FUNBS(p1a+t+(q-1)*(p1a-1),j);
          Sig(q+t,t)+=FUNBS(p1a+t+(q-1)*(p1a-1),j);
          }
      }
    }
  }
  Sig/=k;
  
  /* calculate sigma */
  double epsilon = 0;
  double qq = 0;
  Eigen::MatrixXd bs  = Eigen::MatrixXd::Zero(p1a, p1a);
  Eigen::MatrixXd bs2  = Eigen::MatrixXd::Zero(p1a, p1a);
  Eigen::MatrixXd ZZT = Eigen::MatrixXd::Zero(p1a, p1a);
  
  for (j=0;j<k;j++) 
  {
    
    q=mdata(j);
    
    for(t=0;t<p1a;t++) bs(t, t) = FUNBS(t, j);
    
    if (p1a>1) {
      for(u=1;u<p1a;u++)
      {
        for(t=0;t<(p1a-u);t++) {
          bs(t,u+t) = FUNBS(p1a+t+(u-1)*(p1a-1),j);
          bs(u+t,t) = FUNBS(p1a+t+(u-1)*(p1a-1),j);
        }
      }
    }
    
    for (i=0;i<q;i++) {
      epsilon = Y(mdataS(j)-1+i) - MultVV(X1.row(mdataS(j)-1+i), beta);
      ZZT = MultVVoutprod(Z.row(mdataS(j)-1+i));
      bs2 = ZZT*bs;
      qq += pow(epsilon, 2) - 2*epsilon*MultVV(Z.row(mdataS(j)-1+i), FUNB.col(j)) + bs2.trace();
    }
  }
  
  sigma = 1/mdata.sum()*qq;
  
  /* calculate H01*/
  double dem=0;
  int risk1_index=a-1;
  
  for (j=0;j<k;j++)
  {

    dem+=FUNEC(0,j)*exp(MultVV(X2.row(j), gamma1));
    
    if (cmprsk(j) == 1)
    {
      
      if (j == k-1)
      {
        H01(risk1_index, 2)=H01(risk1_index, 1)/dem;
        risk1_index--;
      }
      else if (survtime(j+1) != survtime(j))
      {
        H01(risk1_index, 2)=H01(risk1_index, 1)/dem;
        risk1_index--;
      }
      
      else
      {
        for (j=j+1;j<k;j++)
        {
          
          dem+=FUNEC(0,j)*exp(MultVV(X2.row(j), gamma1));
          
          if (j == k-1)
          {
            H01(risk1_index, 2)=H01(risk1_index, 1)/dem;
            risk1_index--;
            break;
          }
          else if (survtime(j+1) != survtime(j))
          {
            H01(risk1_index, 2)=H01(risk1_index, 1)/dem;
            risk1_index--;
            break;
          }
          else continue;
        }
      }
      
    }
    else continue;
  }
  
  /* calculate gamma */
  double scalefH01=0;

  risk1_index=a-1;
  
  int p2=gamma1.size();
  
  Eigen::VectorXd SX = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd X22 = Eigen::VectorXd::Zero(p2);
  Eigen::MatrixXd XX = Eigen::MatrixXd::Zero(p2, p2);
  Eigen::MatrixXd SXX = Eigen::MatrixXd::Zero(p2, p2);
  Eigen::VectorXd SX_new = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SX_inter = Eigen::VectorXd::Zero(p2);
  Eigen::MatrixXd SXX_new = Eigen::MatrixXd::Zero(p2, p2);
  
  for (j=0; j<k; j++)
  {
    
    XX = MultVVoutprod(X2.row(j));
    scalef = FUNEC(0,j)*exp(MultVV(X2.row(j), gamma1));
    XX*=scalef;
    SXX+=XX;
    X22=X2.row(j);
    X22*=scalef;
    SX+=X22;
    
    if (cmprsk(j) == 1)
    {
      if (j == k-1)
      {
        scalefH01 = H01(risk1_index, 2);
        SXX*=scalefH01;
        SXX_new+=SXX;
        SXX/=scalefH01;
        SX*=scalefH01;
        SX_new+=SX;
        SX/=scalefH01;
        risk1_index--;
        
      }
      else if (survtime(j+1) != survtime(j))
      {
        scalefH01 = H01(risk1_index, 2);
        SXX*=scalefH01;
        SXX_new+=SXX;
        SXX/=scalefH01;
        SX*=scalefH01;
        SX_new+=SX;
        SX/=scalefH01;
        risk1_index--;
      }
      else
      {
        for (j=j+1;j<k;j++)
        {
          XX = MultVVoutprod(X2.row(j));
          scalef = FUNEC(0,j)*exp(MultVV(X2.row(j), gamma1));
          XX*=scalef;
          SXX+=XX;
          X22=X2.row(j);
          X22*=scalef;
          SX+=X22;
          
          if (j == k-1)
          {
            scalefH01 = H01(risk1_index, 2);
            SXX*=scalefH01;
            SXX_new+=SXX;
            SXX/=scalefH01;
            SX*=scalefH01;
            SX_new+=SX;
            SX/=scalefH01;
            risk1_index--;
            break;
          }
          else if (survtime(j+1) != survtime(j))
          {
            scalefH01 = H01(risk1_index, 2);
            SXX*=scalefH01;
            SXX_new+=SXX;
            SXX/=scalefH01;
            SX*=scalefH01;
            SX_new+=SX;
            SX/=scalefH01;
            risk1_index--;
            break;
          }
          else continue;
        }
      }
    }
    else continue;
  }
  
  for (j=0;j<k;j++)
  {
    if (cmprsk(j) == 1) SX_inter+=X2.row(j);
  }
  
  gamma1+=SXX_new.inverse()*(SX_inter - SX_new);
  
  /*  calculate alpha*/
  
  Eigen::MatrixXd TD = Eigen::MatrixXd::Zero(p1a, p1a);
  Eigen::MatrixXd TDD = Eigen::MatrixXd::Zero(p1a, p1a);
  Eigen::VectorXd TN = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TNN = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd N = Eigen::VectorXd::Zero(p1a);
  
  risk1_index = a-1;
  
  for (j=0;j<k;j++)
  {
    
    for(t=0;t<p1a;t++) bbT(t,t) = FUNBSEC(t,j);
    
    if(p1a>1)
    {
      for(i=1;i<p1a;i++)
      {
        for(t=0;t<p1a-i;t++)   {
          bbT(t,i+t) = FUNBSEC(p1a+t+(i-1)*(p1a-1),j);
          bbT(i+t,t) = bbT(t,i+t);
          }
      }
    }
    
    for (t=0;t<p1a;t++) N(t) = FUNBEC(t,j);
    

    bbT*=exp(MultVV(X2.row(j), gamma1));
    TD+=bbT;
    N*=exp(MultVV(X2.row(j), gamma1));
    TN+=N;
    
    if (cmprsk(j) == 1)
    {
      if (j == k-1)
      {
        TD*=H01(risk1_index,2);
        TDD+=TD;
        TD/=H01(risk1_index,2);
        TN*=H01(risk1_index,2);
        TNN+=TN;
        TN/=H01(risk1_index,2);
        risk1_index--;
        
      }
      else if (survtime(j+1) != survtime(j))
      {
        TD*=H01(risk1_index,2);
        TDD+=TD;
        TD/=H01(risk1_index,2);
        TN*=H01(risk1_index,2);
        TNN+=TN;
        TN/=H01(risk1_index,2);
        risk1_index--;
      }
      else
      {
        for (j=j+1;j<k;j++)
        {
          for(t=0;t<p1a;t++) bbT(t,t) = FUNBSEC(t,j);
          
          if(p1a>1)
          {
            for(i=1;i<p1a;i++)
            {
              for(t=0;t<p1a-i;t++)   {
                bbT(t,i+t) = FUNBSEC(p1a+t+(i-1)*(p1a-1),j);
                bbT(i+t,t) = bbT(t,i+t);
              }
            }
          }
          
          for (t=0;t<p1a;t++) N(t) = FUNBEC(t,j);
          
          bbT*=exp(MultVV(X2.row(j), gamma1));
          TD+=bbT;
          N*=exp(MultVV(X2.row(j), gamma1));
          TN+=N;
          
          if (j == k-1)
          {
            TD*=H01(risk1_index,2);
            TDD+=TD;
            TD/=H01(risk1_index,2);
            TN*=H01(risk1_index,2);
            TNN+=TN;
            TN/=H01(risk1_index,2);
            risk1_index--;
            break;
          }
          else if (survtime(j+1) != survtime(j))
          {
            TD*=H01(risk1_index,2);
            TDD+=TD;
            TD/=H01(risk1_index,2);
            TN*=H01(risk1_index,2);
            TNN+=TN;
            TN/=H01(risk1_index,2);
            risk1_index--;
            break;
          }
          else continue;
        }
      }
    }
    else continue;
  }
  
  TN = Eigen::VectorXd::Zero(p1a);
  
  for (j=0;j<k;j++)
  {
    if(cmprsk(j)==1)
    {
      N = FUNB.col(j);
      TN+=N;
    } else continue;
  }
  
  N=TDD.inverse()*(TN-TNN);
  alpha1+=N;
  
  return Rcpp::List::create(Rcpp::Named("beta")=beta,
                            Rcpp::Named("gamma1")=gamma1,
                            Rcpp::Named("alpha1")=alpha1,
                            Rcpp::Named("H01")= H01,
                            Rcpp::Named("Sig")=Sig,
                            Rcpp::Named("sigma")=sigma);
  
  
  }
