#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Eigen::VectorXd getPus(const Eigen::VectorXd & beta, 
                       const Eigen::VectorXd & gamma1,  
                       const Eigen::VectorXd & alpha1,  
                       const Eigen::MatrixXd & Sig, const double sigma, 
                       const Eigen::MatrixXd & H01,
                       const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
                       const Eigen::VectorXd & Y, const Eigen::MatrixXd & X2,
                       const Eigen::VectorXd & mdata, const Eigen::VectorXd & mdataS,
                       const Eigen::VectorXd & survtime,
                       const Eigen::MatrixXd & xsmatrix, 
                       const Eigen::MatrixXd & wsmatrix,
                       const double u) {
  
  //calculate the square root of random effect covariance matrix 
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd eigenSQ = svd.singularValues();
  int i,j,q,db;
  for (i=0;i<eigenSQ.size();i++) {
    eigenSQ(i) = sqrt(eigenSQ(i));
  }
  Eigen::MatrixXd SigSQRT  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  int p1a=Z.cols();
  int k=X2.rows();
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd weightbi(p1a);
  Eigen::VectorXd ri(p1a);
  Eigen::VectorXd Pus(k);
  
  int point=wsmatrix.rows();
  
  double dems,demu,cuh01,xgamma1,tempu,temps,mu,zb,temp;
  
  double cuh01u=CH(H01, u);
  for (j=0;j<k;j++) {
    dems=0;
    demu=0;
    q=mdata(j);
    cuh01=CH(H01, survtime(j));
    xgamma1=MultVV(X2.row(j),gamma1);
    
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
      
      for (i=0;i<p1a;i++) temp*=weightbi(i);
      
      temps=temp*exp(0-cuh01*exp(xgamma1+MultVV(alpha1,bi)));
      tempu=temp*exp(0-cuh01u*exp(xgamma1+MultVV(alpha1,bi)));
      
      dems+=temps;
      demu+=tempu;
    }
    Pus(j) = demu/dems;
    
  }
  return(Pus);
  
  
  
  }
