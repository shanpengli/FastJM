#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
double getloglikeCstandard(const Eigen::VectorXd & beta, 
                               const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                               const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                               const Eigen::MatrixXd & Sig, const double sigma, 
                               const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
                               const Eigen::VectorXd & Y, const Eigen::MatrixXd & X2, 
                               const Eigen::VectorXd & survtime, 
                               const Eigen::VectorXd & cmprsk, const Eigen::VectorXd & mdata, 
                               const Eigen::VectorXd & mdataS, const Eigen::MatrixXd & xsmatrix, 
                               const Eigen::MatrixXd & wsmatrix, 
                               const Eigen::VectorXd & CUH01, 
                               const Eigen::VectorXd & CUH02, 
                               const Eigen::VectorXd & HAZ01, 
                               const Eigen::VectorXd & HAZ02) {
  
  //calculate the square root of random effect covariance matrix 
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd eigenSQ = svd.singularValues();
  int i,j,q,t,db;
  for (i=0;i<eigenSQ.size();i++) {
    eigenSQ(i) = sqrt(eigenSQ(i));
  }
  Eigen::MatrixXd SigSQRT  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  int k=mdata.size();
  int p1a=Z.cols();
  
  double dem,cuh01,cuh02,haz01,haz02,xgamma1,xgamma2,temp,mu,zb,loglike;
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd weightbi(p1a);
  Eigen::VectorXd ri(p1a);
  
  int point=wsmatrix.rows();
  loglike = 0;
  
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
      temp=1;
      
      for (i=0;i<p1a;i++) bi(i)=ri(i);
      
      for (i=0;i<q;i++) {
        mu=MultVV(X1.row(mdataS(j)-1+i),beta);
        zb=MultVV(Z.row(mdataS(j)-1+i),bi);
        temp*=1/sqrt(sigma*2*M_PI)*exp(-1/(2*sigma)*pow((Y(mdataS(j)-1+i) - mu - zb), 2)); 
      }
      
      if(cmprsk(j)==1)  temp*=haz01*exp(xgamma1+MultVV(alpha1,bi));
      if(cmprsk(j)==2)  temp*=haz02*exp(xgamma2+MultVV(alpha2,bi));
      
      temp*=exp(0-cuh01*exp(xgamma1+MultVV(alpha1,bi))-cuh02*exp(xgamma2+MultVV(alpha2,bi)));
      for (i=0;i<p1a;i++) temp*=weightbi(i);
      temp*=1/sqrt(pow(M_PI, p1a));
      
      dem+=temp;
    }
    loglike+=log(dem);
  }
  
  return loglike;
  
  
  
  
  
  
}
