#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
double getES(const Eigen::VectorXd & beta, const double sigma, 
             const Eigen::VectorXd & gamma1, const Eigen::VectorXd & alpha1, 
             const Eigen::MatrixXd & Sig, 
             const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
             const Eigen::VectorXd & Y, const Eigen::VectorXd & X2, 
             const Eigen::MatrixXd & xsmatrix, 
             const Eigen::MatrixXd & wsmatrix,
             const double CH0s, const double CH0u, const Eigen::VectorXd & Posbi,
             const Eigen::MatrixXd & Poscov){ 
  
  //calculate the square root of random effect covariance matrix 
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig.inverse(), Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd eigenSQ = svd.singularValues();
  int i,j,q,t,db,uu;
  for (i=0;i<eigenSQ.size();i++) {
    eigenSQ(i) = sqrt(eigenSQ(i));
  }
  Eigen::MatrixXd SigSQRT  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  int p1a = Z.cols();
  int ni = Z.rows();
  double xgamma1,temp,mu,zb;
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd bi2(p1a);
  Eigen::VectorXd weightbi(p1a);
  Eigen::VectorXd ri(p1a);
  
  Eigen::VectorXd rii(p1a);
  Eigen::MatrixXd Hi(p1a, p1a);
  
  int point=wsmatrix.rows();
  
  
  xgamma1=MultVV(X2,gamma1);
  
  double S=0;
  double dem=0;
  for (db=0;db<point;db++) {
    //calculate the square root of covariance of Empirical Bayes estimates
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Poscov, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd eigenSQ = svd.singularValues();
    for (i=0;i<eigenSQ.size();i++) {
      eigenSQ(i) = sqrt(eigenSQ(i));
    }
    Hi  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
    
    bi = xsmatrix.row(db);
    weightbi = wsmatrix.row(db);
    ri = Posbi + sqrt(2)*Hi*bi;
    rii = SigSQRT*ri;
    temp=exp(10);
    
    for (i=0;i<p1a;i++) bi(i)=ri(i);
    
    for (i=0;i<ni;i++) {
      mu=MultVV(X1.row(i),beta);
      zb=MultVV(Z.row(i),bi);
      temp*=exp(-1/(2*sigma)*pow((Y(i) - mu - zb), 2));
    }
    
    temp*=exp(0-CH0s*exp(xgamma1+MultVV(alpha1,bi)));
    for (i=0;i<p1a;i++) temp*=weightbi(i);
    bi2 = xsmatrix.row(db);
    temp*=exp(-pow(rii.norm(), 2)/2)*exp(pow(bi2.norm(), 2));
    
    dem+=temp;
    
    //calculate survival
    S+=exp(log(temp) + CH0s*exp(xgamma1+MultVV(alpha1,bi)) - CH0u*exp(xgamma1+MultVV(alpha1,bi)));
    
  }
  
  if(dem==0) {
    Rprintf("Program stops because of the data issue.\n");
    return ( 100.0 );
  }
  
  S/=dem;
  
  return S;
}