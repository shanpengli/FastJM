#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List getECIF(const Eigen::VectorXd & beta, const double sigma, 
                   const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                   const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                   const Eigen::MatrixXd & Sig, 
                   const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
                   const Eigen::VectorXd & Y, const Eigen::VectorXd & X2, 
                   const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                   const Eigen::MatrixXd & xsmatrix, 
                   const Eigen::MatrixXd & wsmatrix,
                   const double CH01, const double CH02,
                   const double s, const double u, const Eigen::VectorXd & Posbi,
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
  
  double CIF1=0;
  double CIF2=0;
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
    
    temp*=exp(0-CH01*exp(xgamma1+MultVV(alpha1,bi))-
      CH02*exp(xgamma2+MultVV(alpha2,bi)));
    for (i=0;i<p1a;i++) temp*=weightbi(i);
    bi2 = xsmatrix.row(db);
    temp*=exp(-pow(rii.norm(), 2)/2)*exp(pow(bi2.norm(), 2));
    
    dem+=temp;
    
    //calculate CIF
    CIF1+=exp(log(temp) + log(GetCIF1CR(gamma1,gamma2,alpha1,alpha2,X2,H01,H02,s,u,ri)) -
      (0-CH01*exp(xgamma1+MultVV(alpha1,bi))-CH02*exp(xgamma2+MultVV(alpha2,bi))));
    CIF2+=exp(log(temp) + log(GetCIF2CR(gamma1,gamma2,alpha1,alpha2,X2,H01,H02,s,u,ri)) -
      (0-CH01*exp(xgamma1+MultVV(alpha1,bi))-CH02*exp(xgamma2+MultVV(alpha2,bi))));
    
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

// [[Rcpp::export]]
Rcpp::List getECIF_JMH(const Eigen::VectorXd & beta, const Eigen::VectorXd & tau, 
                       const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                       const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                       const double nu1, const double nu2, const Eigen::MatrixXd & Sig, 
                       const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
                       const Eigen::MatrixXd & W, const Eigen::VectorXd & Y, 
                       const Eigen::VectorXd & X2, 
                       const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                       const Eigen::MatrixXd & xsmatrix, 
                       const Eigen::MatrixXd & wsmatrix,
                       const double CH01, const double CH02,
                       const double s, const double u){ 
  
  //calculate the square root of random effect covariance matrix 
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd eigenSQ = svd.singularValues();
  int i,j,q,t,db,uu;
  for (i=0;i<eigenSQ.size();i++) {
    eigenSQ(i) = sqrt(eigenSQ(i));
  }
  Eigen::MatrixXd SigSQRT  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  int p1a=Z.cols();
  int ni = Z.rows();
  double xgamma1,xgamma2,temp,mu,sigma,zb,wi;
  Eigen::VectorXd bwi(p1a+1);
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd weightbwi(p1a+1);
  Eigen::VectorXd ri(p1a+1);
  
  int point=wsmatrix.rows();
  
  
  xgamma1=MultVV(X2,gamma1);
  xgamma2=MultVV(X2,gamma2);
  
  double CIF1=0;
  double CIF2=0;
  double dem=0;
  for (db=0;db<point;db++) {
    bwi = xsmatrix.row(db);
    weightbwi = wsmatrix.row(db);
    ri = sqrt(2)*SigSQRT*bwi;
    temp=exp(10);
    
    for (i=0;i<p1a;i++) bi(i)=ri(i);
    wi=ri(p1a);
    
    for (i=0;i<ni;i++) {
      mu=MultVV(X1.row(i),beta);
      zb=MultVV(Z.row(i),bi);
      sigma=exp(MultVV(W.row(i),tau) + wi);
      temp*=1/sqrt(sigma)*exp(-1/(2*sigma)*pow((Y(i) - mu - zb), 2));
    }
    
    temp*=exp(0-CH01*exp(xgamma1+MultVV(alpha1,bi)+nu1*wi)-
      CH02*exp(xgamma2+MultVV(alpha2,bi)+nu2*wi));
    for (i=0;i<(p1a+1);i++) temp*=weightbwi(i);
    
    dem+=temp;
    
    //calculate CIF
    CIF1+=exp(log(temp) + log(GetCIF1CR_JMH(gamma1,gamma2,alpha1,alpha2,nu1,nu2,X2,H01,H02,s,u,ri,(p1a+1))) -
      (0-CH01*exp(xgamma1+MultVV(alpha1,bi)+nu1*wi)-CH02*exp(xgamma2+MultVV(alpha2,bi)+nu2*wi)));
    CIF2+=exp(log(temp) + log(GetCIF2CR_JMH(gamma1,gamma2,alpha1,alpha2,nu1,nu2,X2,H01,H02,s,u,ri,(p1a+1))) -
      (0-CH01*exp(xgamma1+MultVV(alpha1,bi)+nu1*wi)-CH02*exp(xgamma2+MultVV(alpha2,bi)+nu2*wi)));
    
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

// [[Rcpp::export]]
Rcpp::List getECIFad_JMH(const Eigen::VectorXd & beta, const Eigen::VectorXd & tau, 
                         const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                         const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                         const double nu1, const double nu2, const Eigen::MatrixXd & Sig, 
                         const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
                         const Eigen::MatrixXd & W, const Eigen::VectorXd & Y, 
                         const Eigen::VectorXd & X2, 
                         const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                         const Eigen::MatrixXd & xsmatrix, 
                         const Eigen::MatrixXd & wsmatrix,
                         const double CH01, const double CH02,
                         const double s, const double u,
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
  
  int p1a=Z.cols();
  int ni = Z.rows();
  double xgamma1,xgamma2,temp,mu,sigma,zb,wi;
  Eigen::VectorXd bwi(p1a+1);
  Eigen::VectorXd bwii(p1a+1);
  Eigen::VectorXd bwi2(p1a+1);
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd weightbwi(p1a+1);
  Eigen::VectorXd ri(p1a+1);
  Eigen::VectorXd rii(p1a+1);
  
  Eigen::MatrixXd Hi(p1a+1, p1a+1);
  Eigen::MatrixXd Hi2(p1a+1, p1a+1);
  
  int point=wsmatrix.rows();
  
  
  xgamma1=MultVV(X2,gamma1);
  xgamma2=MultVV(X2,gamma2);
  
  double CIF1=0;
  double CIF2=0;
  double dem=0;
  
  Hi = Poscov;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd2(Hi, Eigen::ComputeThinU | Eigen::ComputeThinV);
  eigenSQ = svd2.singularValues();
  for (i=0;i<eigenSQ.size();i++) {
    eigenSQ(i) = sqrt(eigenSQ(i));
  }
  Hi2  = svd2.matrixU() * eigenSQ.asDiagonal() * svd2.matrixV().transpose();
  
  
  for (db=0;db<point;db++) {
    bwi = xsmatrix.row(db);
    weightbwi = wsmatrix.row(db);
    ri = Posmean + sqrt(2)*Hi2*bwi;
    rii = SigSQRT*ri;
    temp=exp(10);
    
    for (i=0;i<p1a;i++) bi(i)=ri(i);
    wi=ri(p1a);
    
    for (i=0;i<ni;i++) {
      mu=MultVV(X1.row(i),beta);
      zb=MultVV(Z.row(i),bi);
      sigma=exp(MultVV(W.row(i),tau) + wi);
      temp*=1/sqrt(sigma)*exp(-1/(2*sigma)*pow((Y(i) - mu - zb), 2));
    }
    
    temp*=exp(0-CH01*exp(xgamma1+MultVV(alpha1,bi)+nu1*wi)-
      CH02*exp(xgamma2+MultVV(alpha2,bi)+nu2*wi));
    for (i=0;i<(p1a+1);i++) temp*=weightbwi(i);
    bwi2 = xsmatrix.row(db);
    temp*=exp(-pow(rii.norm(), 2)/2)*exp(pow(bwi2.norm(), 2));
    
    dem+=temp;
    
    //calculate CIF
    CIF1+=exp(log(temp) + log(GetCIF1CR_JMH(gamma1,gamma2,alpha1,alpha2,nu1,nu2,X2,H01,H02,s,u,ri,(p1a+1))) -
      (0-CH01*exp(xgamma1+MultVV(alpha1,bi)+nu1*wi)-CH02*exp(xgamma2+MultVV(alpha2,bi)+nu2*wi)));
    CIF2+=exp(log(temp) + log(GetCIF2CR_JMH(gamma1,gamma2,alpha1,alpha2,nu1,nu2,X2,H01,H02,s,u,ri,(p1a+1))) -
      (0-CH01*exp(xgamma1+MultVV(alpha1,bi)+nu1*wi)-CH02*exp(xgamma2+MultVV(alpha2,bi)+nu2*wi)));
    
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
