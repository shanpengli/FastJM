#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List getCovSF(const Eigen::VectorXd & beta, const Eigen::VectorXd & gamma1, 
                  const Eigen::VectorXd & alpha1, 
                  const Eigen::MatrixXd & H01, 
                  const Eigen::MatrixXd & Sig, 
                  const double sigma, const Eigen::MatrixXd & Z, 
                  const Eigen::MatrixXd & X1, const Eigen::VectorXd & Y, 
                  const Eigen::MatrixXd & X2, const Eigen::VectorXd & survtime, 
                  const Eigen::VectorXd & cmprsk, const Eigen::VectorXd & mdata, 
                  const Eigen::VectorXd & mdataS,
                  const Eigen::MatrixXd & FUNBS, 
                  const Eigen::MatrixXd & FUNEC, 
                  const Eigen::MatrixXd & FUNBEC,
                  const Eigen::MatrixXd & FUNBSEC, 
                  const Eigen::MatrixXd & FUNB) {
  
  int d = beta.size() + gamma1.size() + alpha1.size() + Sig.cols()*(Sig.cols() + 1)/2 + 1;
  int a = H01.rows();
  int k = mdata.size();
  
  int i,q,j,t;
  
  double temp,temp1;
  
  Eigen::MatrixXd SS  = Eigen::MatrixXd::Zero(d,d);
  Eigen::MatrixXd SSinv  = Eigen::MatrixXd::Zero(d,d);
  Eigen::VectorXd S = Eigen::VectorXd::Zero(d);
  
  int risk1_index;
  int risk1_index_temp=a-1;
  int risk1_index_ttemp=a-1;
  int risk1_index_tttemp=a-1;
  int risk1_index_vtemp=a-1;
  int risk1_index_vttemp=a-1;
  int risk1_index_vtttemp=a-1;
  
  Eigen::VectorXd CumuH01 = Eigen::VectorXd::Zero(a);
  
  temp1=0;
  for (j=0;j<a;j++) {
    temp1 += H01(j, 2);
    CumuH01(j) = temp1;
    }
  
  double epsilon=0;
  double qqsigma=0;
  int p1 = beta.size();
  int p1a = Z.cols();
  int p2 = gamma1.size();
  
  Eigen::VectorXd SZ = Eigen::VectorXd::Zero(p1);
  Eigen::VectorXd SZ1 = Eigen::VectorXd::Zero(p1);
  Eigen::MatrixXd SZZ = Eigen::MatrixXd::Zero(p1, p1a);
  
  Eigen::VectorXd X = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SX = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SX1 = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SRX = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SRXX = Eigen::VectorXd::Zero(p2);
  Eigen::MatrixXd SXX1 = Eigen::MatrixXd::Zero(p2, a);
  Eigen::MatrixXd SXX11 = Eigen::MatrixXd::Zero(p2, a);
  Eigen::MatrixXd SRXX1 = Eigen::MatrixXd::Zero(p2, k);
  
  Eigen::VectorXd N = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TN = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TN1 = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TRN = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TRNN = Eigen::VectorXd::Zero(p1a);
  Eigen::MatrixXd TRNN1 = Eigen::MatrixXd::Zero(p1a,k);
  Eigen::MatrixXd TNN1 = Eigen::MatrixXd::Zero(p1a,a);
  Eigen::MatrixXd TNN11 = Eigen::MatrixXd::Zero(p1a,a);
  
  Eigen::MatrixXd bs = Eigen::MatrixXd::Zero(p1a,p1a);
  Eigen::MatrixXd ZZT = Eigen::MatrixXd::Zero(p1a,p1a);
  Eigen::MatrixXd bs2 = Eigen::MatrixXd::Zero(p1a,p1a);
  
  for (j=0;j<k;j++) {
    
    S = Eigen::VectorXd::Zero(d);
    q = mdata(j);
    
    /* calculate score for beta */
    SZ = Eigen::VectorXd::Zero(p1);
    for (i=0;i<q;i++) {
      
      epsilon = Y(mdataS(j)-1+i) - MultVV(X1.row(mdataS(j)-1+i), beta) - 
        MultVV(Z.row(mdataS(j)-1+i), FUNB.col(j));
      SZ += epsilon*X1.row(mdataS(j)-1+i);
      
    }
    SZ /= sigma;
    
    for (i=0;i<p1;i++) S(i) = SZ(i);
    
    /* calculate score for gamma */
    if (j == 0)
    {
      temp=0;
      risk1_index=risk1_index_temp;
      for (q=j;q<k;q++)
      {
        
        temp+=exp(MultVV(X2.row(q), gamma1))*FUNEC(0,q);
        SX += exp(MultVV(X2.row(q), gamma1))*FUNEC(0,q)*X2.row(q);
        
        if (cmprsk(q) == 1)
        {
          if (q == k-1)
          {
            SX*= H01(risk1_index, 1)/pow(temp, 2);
            SXX1.col(a-1-risk1_index) = SX;
            SRXX += SX;
            SX/= H01(risk1_index, 1)/pow(temp, 2);
            SX/=temp;
            SXX11.col(a-1-risk1_index) = SX;
            SX*=temp;
            risk1_index--;
          }
          else if (survtime(q+1) != survtime(q))
          {
            SX*= H01(risk1_index, 1)/pow(temp, 2);
            SXX1.col(a-1-risk1_index) = SX;
            SRXX += SX;
            SX/= H01(risk1_index, 1)/pow(temp, 2);
            SX/=temp;
            SXX11.col(a-1-risk1_index) = SX;
            SX*=temp;
            risk1_index--;
          }
          else
          {
            for (q=q+1;q<k;q++)
            {
              temp+=exp(MultVV(X2.row(q), gamma1))*FUNEC(0,q);
              SX += exp(MultVV(X2.row(q), gamma1))*FUNEC(0,q)*X2.row(q);
              
              if (q == k-1)
              {
                SX*= H01(risk1_index, 1)/pow(temp, 2);
                SXX1.col(a-1-risk1_index) = SX;
                SRXX += SX;
                SX/= H01(risk1_index, 1)/pow(temp, 2);
                SX/=temp;
                SXX11.col(a-1-risk1_index) = SX;
                SX*=temp;
                risk1_index--;
                break;
              }
              else if (survtime(q+1) != survtime(q))
              {
                SX*= H01(risk1_index, 1)/pow(temp, 2);
                SXX1.col(a-1-risk1_index) = SX;
                SRXX += SX;
                SX/= H01(risk1_index, 1)/pow(temp, 2);
                SX/=temp;
                SXX11.col(a-1-risk1_index) = SX;
                SX*=temp;
                risk1_index--;
                break;
              }
              else continue;
            }
          }
          
        }
        else continue;
      }
      SRXX1.col(j) = SRXX;
    }
    else
    {
      if (risk1_index_temp>=0)
      {
        if (survtime(j) >= H01(risk1_index_temp, 0))
        {
          SRXX1.col(j) = SRXX1.col(j-1);
        }
        else
        {
          risk1_index_temp--;
          if (risk1_index_temp>=0)
          {
            SRXX = SRXX1.col(j-1);
            SRXX -= SXX1.col(a-1-risk1_index_temp-1);
            SRXX1.col(j) = SRXX;
          }
        }
      }
      else
      {
        risk1_index_temp=0;
      }
    }
    SRX = SRXX1.col(j);
    
    if (j==0)
    {
      SRX -= CumuH01(risk1_index_ttemp)*X2.row(j);  
    }
    else if (survtime(j) >= H01(risk1_index_ttemp, 0))
    {
      SRX -= CumuH01(risk1_index_ttemp)*X2.row(j);  
    }
    else
    {
      risk1_index_ttemp--;
      if (risk1_index_ttemp>=0)
      {
        SRX -= CumuH01(risk1_index_ttemp)*X2.row(j);  
      }
      else
      {
        SRX = SRXX1.col(j);
        risk1_index_ttemp=0;
      }
    }
    
    SRX*= exp(MultVV(X2.row(j), gamma1))*FUNEC(0,j);
    
    if (survtime(j) >= H01(risk1_index_tttemp, 0))
    {
      if (cmprsk(j) == 1)
      {
        X = X2.row(j);
        X -= SXX11.col(a-1-risk1_index_tttemp);
        X += SRX;
        for (q=0;q<p2;q++) S(p1+q) = X(q);
      }
      else
      {
        for (q=0;q<p2;q++) S(p1+q) = SRX(q);
      }
    }
    else
    {
      risk1_index_tttemp--;
      if (risk1_index_tttemp>=0)
      {
        if (cmprsk(j) == 1)
        {
          X = X2.row(j);
          X -= SXX11.col(a-1-risk1_index_tttemp);
          X += SRX;
          for (q=0;q<p2;q++) S(p1+q) = X(q);
        }
        else
        {
          for (q=0;q<p2;q++) S(p1+q) = SRX(q);
        }
      }
      else
      {
        risk1_index_tttemp=0;
        for (q=0;q<p2;q++) S(p1+q) = SRX(q);
      }
    }
    
    /* calculate score for alpha */
    /*  alpha1 */
    if (j == 0)
    {
      temp=0;
      
      TN = Eigen::VectorXd::Zero(p1a);
      TRN = Eigen::VectorXd::Zero(p1a);
      
      risk1_index=risk1_index_vtemp;
      for (q=j;q<k;q++)
      {
        temp += exp(MultVV(X2.row(q), gamma1))*FUNEC(0, q);
        for (i=0;i<p1a;i++) N(i) = FUNBEC(i,q);
        TN += exp(MultVV(X2.row(q), gamma1))*N;
        if (cmprsk(q) == 1)
        {
          if (q == k-1)
          {
            TN *= H01(risk1_index, 1)/pow(temp,2);
            TNN1.col(a-1-risk1_index) = TN;
            TRNN += TN;
            TN /= H01(risk1_index, 1)/pow(temp,2);
            TN /= temp;
            TNN11.col(a-1-risk1_index) = TN;
            TN *= temp;
            risk1_index--;
          }
          else if (survtime(q+1) != survtime(q))
          {
            TN *= H01(risk1_index, 1)/pow(temp,2);
            TNN1.col(a-1-risk1_index) = TN;
            TRNN += TN;
            TN /= H01(risk1_index, 1)/pow(temp,2);
            TN /= temp;
            TNN11.col(a-1-risk1_index) = TN;
            TN *= temp;
            risk1_index--;
          }
          else
          {
            for (q=q+1;q<k;q++)
            {
              temp += exp(MultVV(X2.row(q), gamma1))*FUNEC(0, q);
              for (i=0;i<p1a;i++) N(i) = FUNBEC(i,q);
              TN += exp(MultVV(X2.row(q), gamma1))*N;
              if (q == k-1)
              {
                TN *= H01(risk1_index, 1)/pow(temp,2);
                TNN1.col(a-1-risk1_index) = TN;
                TRNN += TN;
                TN /= H01(risk1_index, 1)/pow(temp,2);
                TN /= temp;
                TNN11.col(a-1-risk1_index) = TN;
                TN *= temp;
                risk1_index--;
                break;
              }
              else if (survtime(q+1) != survtime(q))
              {
                TN *= H01(risk1_index, 1)/pow(temp,2);
                TNN1.col(a-1-risk1_index) = TN;
                TRNN += TN;
                TN /= H01(risk1_index, 1)/pow(temp,2);
                TN /= temp;
                TNN11.col(a-1-risk1_index) = TN;
                TN *= temp;
                risk1_index--;
                break;
              }
              else continue;
            }
          }
          
        }
        else continue;
      }
      TRNN1.col(j) = TRNN;
    }
    else
    {
      if (risk1_index_vtemp>=0)
      {
        if (survtime(j) >= H01(risk1_index_vtemp, 0))
        {
          TRNN1.col(j) = TRNN1.col(j-1);
        }
        else
        {
          risk1_index_vtemp--;
          if (risk1_index_vtemp>=0)
          {
            TRNN = TRNN1.col(j-1);
            TRNN -= TNN1.col(a-1-risk1_index_vtemp-1);
            TRNN1.col(j) = TRNN;
          }
        }
      }
      else
      {
        risk1_index_vtemp=0;
      }
    }
    TRN = TRNN1.col(j);
    
    TRN *= exp(MultVV(X2.row(j), gamma1))*FUNEC(0,j);
    
    if (j==0)
    {
      for (t=0;t<p1a;t++) N(t) = FUNBEC(t,j);
      N *= CumuH01(risk1_index_vttemp)*exp(MultVV(X2.row(j), gamma1));
      TRN -= N;
    }
    else if (survtime(j) >= H01(risk1_index_vttemp,0))
    {
      for (t=0;t<p1a;t++) N(t) = FUNBEC(t,j);
      N *= CumuH01(risk1_index_vttemp)*exp(MultVV(X2.row(j), gamma1));
      TRN -= N;
    }
    else
    {
      risk1_index_vttemp--;
      if (risk1_index_vttemp>=0)
      {
        for (t=0;t<p1a;t++) N(t) = FUNBEC(t,j);
        N *= CumuH01(risk1_index_vttemp)*exp(MultVV(X2.row(j), gamma1));
        TRN -= N;
      }
      else
      {
        risk1_index_vttemp=0;
      }
    }
    
    
    if (survtime(j) >= H01(risk1_index_vtttemp,0))
    {
      if (cmprsk(j) == 1)
      {
        TN = FUNB.col(j) - TNN11.col(a-1-risk1_index_vtttemp);
        TN += TRN;
        for (q=0;q<p1a;q++) S(p1+p2+q) = TN(q);
      }
      else
      {
        for (q=0;q<p1a;q++) S(p1+p2+q) = TRN(q);
      }
    }
    else
    {
      risk1_index_vtttemp--;
      if (risk1_index_vtttemp>=0)
      {
        if (cmprsk(j) == 1)
        {
          TN = FUNB.col(j) - TNN11.col(a-1-risk1_index_vtttemp);
          TN += TRN;
          for (q=0;q<p1a;q++) S(p1+p2+q) = TN(q);
        }
        else
        {
          for (q=0;q<p1a;q++) S(p1+p2+q) = TRN(q);
        }
      }
      else
      {
        risk1_index_vtttemp=0;
        for (q=0;q<p1a;q++) S(p1+p2+q) = TRN(q);
      }
    }
    
    /*calculate sigma*/
    for (i=0;i<p1a;i++) bs(i,i) = FUNBS(i, j);
    if(p1a>1)
    {
      for(i=1;i<p1a;i++)
      {
        for(t=0;t<p1a-i;t++) {
          bs(t,i+t) = FUNBS(p1a+t+(i-1)*(p1a-1),j);
          bs(i+t,t) = bs(t,i+t);
        }
      }
    }
    qqsigma = 0;
    q = mdata(j);
    for (i=0;i<q;i++) {
      
      epsilon = Y(mdataS(j)-1+i) - MultVV(X1.row(mdataS(j)-1+i), beta);
      ZZT = MultVVoutprod(Z.row(mdataS(j)-1+i));
      bs2 = ZZT*bs;
      qqsigma += pow(epsilon, 2) + bs2.trace() - 2*epsilon*MultVV(Z.row(mdataS(j)-1+i), FUNB.col(j));
      
    }
    qqsigma /= 2*pow(sigma, 2);
    qqsigma -= q/(2*sigma);
    S(p1+p2+p1a) = qqsigma;
    
    /*calculate Sig matrix*/
    
    for(t=0;t<p1a;t++) bs(t,t) = FUNBS(t,j);
    
    if(p1a>1)
    {
      for(i=1;i<p1a;i++)
      {
        for(t=0;t<p1a-i;t++) {
          bs(t,i+t) = FUNBS(p1a+t+(i-1)*(p1a-1),j);
          bs(i+t,t) = bs(t,i+t);
        }
      }
    }
    
    bs = Sig.inverse()*bs*Sig.inverse() - Sig.inverse();
    
    for (t=0;t<p1a;t++) S(p1+p2+p1a+1+t) = 0.5*bs(t,t);
    
    for(q=1;q<p1a;q++)
    {
      for(t=0;t<(p1a-q);t++) S(p1+p2+p1a+1+p1a+t+(q-1)*(p1a-1)) = bs(t,q+t);
    }
    
    SS += MultVVoutprod(S);
  }
  
  SSinv = SS.inverse();
  
  Eigen::VectorXd sebeta = Eigen::VectorXd::Zero(p1);
  Eigen::VectorXd segamma1 = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd sealpha1 = Eigen::VectorXd::Zero(p1a);
  double sesigma=0;
  Eigen::MatrixXd seSig = Eigen::MatrixXd::Zero(p1a, p1a);
  
  for (t=0;t<p1;t++) sebeta(t) = sqrt(SSinv(t,t));
  for (t=0;t<p2;t++) segamma1(t) = sqrt(SSinv(p1+t,p1+t));
  for (t=0;t<p1a;t++) sealpha1(t) = sqrt(SSinv(p1+p2+t,p1+p2+t));
  sesigma = sqrt(SSinv(p1+p2+p1a,p1+p2+p1a));
  for (t=0;t<p1a;t++) seSig(t,t) = sqrt(SSinv(p1+p2+p1a+1+t,p1+p2+p1a+1+t));
  for(q=1;q<p1a;q++)
  {
    for(t=0;t<(p1a-q);t++) {
      seSig(t,q+t) = sqrt(SSinv(p1+p2+p1a+1+p1a+t+(q-1)*(p1a-1),p1+p2+p1a+1+p1a+t+(q-1)*(p1a-1)));
      seSig(q+t,t) = seSig(t,q+t);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("vcov")=SSinv,
                            Rcpp::Named("sebeta")=sebeta,
                            Rcpp::Named("sesigma")=sesigma,
                            Rcpp::Named("segamma1")=segamma1,
                            Rcpp::Named("sealpha1")=sealpha1,
                            Rcpp::Named("seSig")=seSig);
  
}
