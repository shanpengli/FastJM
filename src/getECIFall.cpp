#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List getECIFall(const Eigen::VectorXd & beta, const double sigma,  
                      const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                      const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                      const Eigen::MatrixXd & Sig, 
                      const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
                      const Eigen::VectorXd & Y, 
                      const Eigen::VectorXd & X2, 
                      const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                      const Eigen::MatrixXd & xsmatrix, 
                      const Eigen::MatrixXd & wsmatrix,
                      const double CH01, const double CH02,
                      const double s, const Eigen::VectorXd & timecif1,
                      const Eigen::VectorXd & timecif2,
                      const Eigen::VectorXd & Posmean,
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
  double xgamma1,xgamma2,temp,mu,zb;
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd bi2(p1a);
  Eigen::VectorXd weightbi(p1a);
  Eigen::VectorXd ri(p1a);
  
  Eigen::VectorXd rii(p1a);
  Eigen::MatrixXd Hi(p1a, p1a);
  
  int point=wsmatrix.rows();
  
  
  xgamma1=MultVV(X2,gamma1);
  xgamma2=MultVV(X2,gamma2);
  
  int a1 = timecif1.size();
  int b1 = timecif2.size();
  
  Eigen::VectorXd CIF1 = Eigen::VectorXd::Zero(a1);
  Eigen::VectorXd CIF2 = Eigen::VectorXd::Zero(b1);
  double dem=0;

  
  for (db=0;db<point;db++) {
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Poscov, Eigen::ComputeThinU | Eigen::ComputeThinV);
    eigenSQ = svd.singularValues();
    for (i=0;i<eigenSQ.size();i++) {
      eigenSQ(i) = sqrt(eigenSQ(i));
    }
    Hi  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
    
    bi = xsmatrix.row(db);
    weightbi = wsmatrix.row(db);
    ri = Posmean + sqrt(2)*Hi*bi;
    rii = SigSQRT*ri;
    temp=exp(10);
    
    for (i=0;i<p1a;i++) bi(i)=ri(i);
    
    for (i=0;i<ni;i++) {
      mu=MultVV(X1.row(i),beta);
      zb=MultVV(Z.row(i),bi);
      temp*=exp(-1/(2*sigma)*pow((Y(i) - mu - zb), 2));
    }
    
    temp*=exp(0-CH01*exp(xgamma1+MultVV(alpha1,bi))-
      CH02*exp(xgamma2+MultVV(alpha2,bi)));
    for (i=0;i<p1a;i++) temp*=weightbi(i);
    bi2 = xsmatrix.row(db);
    temp*=exp(-pow(rii.norm(), 2)/2)*exp(pow(bi2.norm(), 2));
    
    dem+=temp;
    
    Eigen::VectorXd GetCIF1all = GetCIF1CRall(gamma1,gamma2,alpha1,alpha2,X2,H01,H02,s,timecif1,ri);
    Eigen::VectorXd GetCIF2all = GetCIF2CRall(gamma1,gamma2,alpha1,alpha2,X2,H01,H02,s,timecif2,ri);
    
    //calculate CIF across all uncensored event times
    for (i=0;i<a1;i++) {
      CIF1(i)+=exp(log(temp) + log(GetCIF1all(i)) -
        (0-CH01*exp(xgamma1+MultVV(alpha1,bi))-CH02*exp(xgamma2+MultVV(alpha2,bi))));
    }
    
    for (i=0;i<b1;i++) {
      CIF2(i)+=exp(log(temp) + log(GetCIF2all(i)) -
        (0-CH01*exp(xgamma1+MultVV(alpha1,bi))-CH02*exp(xgamma2+MultVV(alpha2,bi))));
    }
    
  }
  
  if(dem==0) {
    Rprintf("Program stops because of the data issue.\n");
    return ( 100.0 );
  }
  
  CIF1/=dem;
  CIF2/=dem;
  
  return Rcpp::List::create(Rcpp::Named("CIF1")=CIF1,
                            Rcpp::Named("CIF2")=CIF2);
}