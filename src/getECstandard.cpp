#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List getECstandard(const Eigen::VectorXd & beta, 
           const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
           const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
           const Eigen::MatrixXd & Sig, const double sigma, 
           const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
           const Eigen::VectorXd & Y, const Eigen::MatrixXd & X2, 
           const Eigen::VectorXd & survtime, 
           const Eigen::VectorXd & cmprsk, const Eigen::VectorXd & mdata, 
           const Eigen::VectorXi & mdataS, const Eigen::MatrixXd & xsmatrix, 
           const Eigen::MatrixXd & wsmatrix, 
           const Eigen::VectorXd & CUH01, 
           const Eigen::VectorXd & CUH02, 
           const Eigen::VectorXd & HAZ01, 
           const Eigen::VectorXd & HAZ02){ 
  
  //calculate the square root of random effect covariance matrix 
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd eigenSQ = svd.singularValues();
  int i,j,q,t,db,u;
  for (i=0;i<eigenSQ.size();i++) {
    eigenSQ(i) = sqrt(eigenSQ(i));
    }
  Eigen::MatrixXd SigSQRT  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  int k=mdata.size();
  int p1a=Z.cols();
  
  double dem,cuh01,cuh02,haz01,haz02,xgamma1,xgamma2,temp,mu,zb;
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd weightbi(p1a);
  Eigen::VectorXd ri(p1a);
  
  /* Define functions of bi*/
  //b
  Eigen::MatrixXd FUNB = Eigen::MatrixXd::Zero(p1a,k);
  //bbT
  Eigen::MatrixXd FUNBS = Eigen::MatrixXd::Zero(p1a*(p1a+1)/2,k);
  //exp(alpha*b)
  Eigen::MatrixXd FUNEC = Eigen::MatrixXd::Zero(2,k);
  //bexp(alpha*b)
  Eigen::MatrixXd FUNBEC = Eigen::MatrixXd::Zero(2*p1a,k);
  //bbTexp(alpha*b)
  Eigen::MatrixXd FUNBSEC = Eigen::MatrixXd::Zero(2*p1a*(p1a+1)/2,k);
  
  int point=wsmatrix.rows();
  
  for(j=0;j<k;j++)
  {
    dem=0;
    q=mdata(j);
    cuh01=CUH01(j);
    cuh02=CUH02(j);
    haz01=HAZ01(j);
    haz02=HAZ02(j);
  
    xgamma1=MultVV(X2.row(j),gamma1);
    
    xgamma2=MultVV(X2.row(j),gamma2);
    
    for (db=0;db<point;db++) {
      
      bi = xsmatrix.row(db);
      weightbi = wsmatrix.row(db);
      ri = sqrt(2)*SigSQRT*bi;
      temp=exp(10);
      
      for (i=0;i<p1a;i++) bi(i)=ri(i);
      
      for (i=0;i<q;i++) {
        mu=MultVV(X1.row(mdataS(j)-1+i),beta);
        zb=MultVV(Z.row(mdataS(j)-1+i),bi);
        temp*=exp(-1/(2*sigma)*pow((Y(mdataS(j)-1+i) - mu - zb), 2)); 
        }
      
      if(cmprsk(j)==1)  temp*=haz01*exp(xgamma1+MultVV(alpha1,bi));
      if(cmprsk(j)==2)  temp*=haz02*exp(xgamma2+MultVV(alpha2,bi));
      
      temp*=exp(0-cuh01*exp(xgamma1+MultVV(alpha1,bi))-cuh02*exp(xgamma2+MultVV(alpha2,bi)));
      for (i=0;i<p1a;i++) temp*=weightbi(i);
      
      dem+=temp;
      
      //calculate h(bi)
      FUNB.col(j)+=temp*bi;
      
      for (i=0;i<p1a;i++) {
        FUNBS(i,j)+=temp*pow(bi(i),2);
        }
      
      if (p1a >= 2) {
        u=0;
        for(i=1;i<p1a;i++)
        {
          for(t=0;t<p1a-i;t++) {
            FUNBS(p1a+u,j) += temp*bi(t)*bi(t+i);
            u++;
          }   
        }
        
      }
      
      FUNEC(0,j)+=temp*exp(MultVV(alpha1,bi));
      FUNEC(1,j)+=temp*exp(MultVV(alpha2,bi));
      
      for (i=0;i<p1a;i++) {
        FUNBEC(i,j)+=temp*bi(i)*exp(MultVV(alpha1,bi));
        FUNBEC(p1a+i,j)+=temp*bi(i)*exp(MultVV(alpha2,bi));
        }
      
      for (i=0;i<p1a;i++) {
        FUNBSEC(i,j)+=temp*exp(MultVV(alpha1,bi))*pow(bi(i),2);
        FUNBSEC(p1a*(p1a+1)/2+i,j)+=temp*exp(MultVV(alpha2,bi))*pow(bi(i),2);
        }
      
      if (p1a >= 2) {
        u=0;
        for(i=1;i<p1a;i++)
        {
          for(t=0;t<p1a-i;t++)
          {
            FUNBSEC(p1a+u,j)+=temp*exp(MultVV(alpha1,bi))*bi(t)*bi(t+i);
            FUNBSEC(p1a*(p1a+1)/2+p1a+u,j)+=temp*exp(MultVV(alpha2,bi))*bi(t)*bi(t+i);
            u++;
          }
        }
      }
      
      
      }
      
    if(dem==0) {
      Rprintf("E step ran into issue for the %dth subject. Program stops.\n", j);
      return ( 100.0 );
    } 
    
    FUNB.col(j)/=dem;
    FUNBS.col(j)/=dem;
    FUNEC.col(j)/=dem;
    FUNBEC.col(j)/=dem;
    FUNBSEC.col(j)/=dem;

  }
  
  return Rcpp::List::create(Rcpp::Named("FUNB")=FUNB,
                            Rcpp::Named("FUNBS")=FUNBS,
                            Rcpp::Named("FUNEC")=FUNEC,
                            Rcpp::Named("FUNBEC")=FUNBEC,
                            Rcpp::Named("FUNBSEC")=FUNBSEC);
}


