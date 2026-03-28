#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
double getloglikeCSFad(const Eigen::VectorXd & beta, const Eigen::VectorXd & tau, 
                       const Eigen::VectorXd & gamma1,
                       const Eigen::VectorXd & alpha1,
                       const double vee1, const Eigen::MatrixXd & Sig, 
                       const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
                       const Eigen::MatrixXd & W, const Eigen::VectorXd & Y, 
                       const Eigen::MatrixXd & X2, const Eigen::VectorXd & survtime, 
                       const Eigen::VectorXd & cmprsk, const Eigen::VectorXd & mdata, 
                       const Eigen::VectorXi & mdataS, const Eigen::MatrixXd & xsmatrix, 
                       const Eigen::MatrixXd & wsmatrix, 
                       const Eigen::VectorXd & CUH01, 
                       const Eigen::VectorXd & HAZ01, 
                       const Eigen::MatrixXd & Posbwi,
                       const Eigen::MatrixXd & Poscov){ 
  
  //calculate the square root of random effect covariance matrix 
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig.inverse(), Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd eigenSQ = svd.singularValues();
  int i,j,q,t,db;
  for (i=0;i<eigenSQ.size();i++) {
    eigenSQ(i) = sqrt(eigenSQ(i));
  }
  Eigen::MatrixXd SigSQRT  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  int k=mdata.size();
  int p1a=Z.cols();
  
  double dem,cuh01,haz01,xgamma1,temp,mu,sigma,zb,wi;
  Eigen::VectorXd bwi(p1a+1);
  Eigen::VectorXd bwii(p1a+1);
  Eigen::VectorXd bwi2(p1a+1);
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd weightbwi(p1a+1);
  Eigen::VectorXd ri(p1a+1);
  Eigen::VectorXd rii(p1a+1);
  
  Eigen::MatrixXd Hi(p1a+1, p1a+1);
  Eigen::MatrixXd Hi2(p1a+1, p1a+1);
  double loglike = 0;
  
  int point=wsmatrix.rows();
  
  for(j=0;j<k;j++)
  {
    dem=0;
    q=mdata(j);
    cuh01=CUH01(j);
    haz01=HAZ01(j);
    
    xgamma1=MultVV(X2.row(j),gamma1);
    
    //calculate the square root of covariance of Empirical Bayes estimates
    for (i=0;i<(p1a+1);i++) Hi.row(i) = Poscov.row(j*(p1a+1)+i);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Hi, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd eigenSQ = svd.singularValues();
    for (i=0;i<eigenSQ.size();i++) {
      eigenSQ(i) = sqrt(eigenSQ(i));
    }
    Hi2  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
    
    for (db=0;db<point;db++) {
      
      bwi = xsmatrix.row(db);
      weightbwi = wsmatrix.row(db);
      bwii = Posbwi.row(j);
      ri = bwii + sqrt(2)*Hi2*bwi;
      rii = SigSQRT*ri;
      temp=1;
      
      for (i=0;i<p1a;i++) bi(i)=ri(i);
      wi=ri(p1a);
      
      for (i=0;i<q;i++) {
        mu=MultVV(X1.row(mdataS(j)-1+i),beta);
        zb=MultVV(Z.row(mdataS(j)-1+i),bi);
        sigma=exp(MultVV(W.row(mdataS(j)-1+i),tau) + wi);
        temp*=1/sqrt(sigma*2*M_PI)*exp(-1/(2*sigma)*pow((Y(mdataS(j)-1+i) - mu - zb), 2)); 
      }
      
      if(cmprsk(j)==1)  temp*=haz01*exp(xgamma1+MultVV(alpha1,bi)+vee1*wi);
      
      temp*=exp(0-cuh01*exp(xgamma1+MultVV(alpha1,bi)+vee1*wi));
      for (i=0;i<(p1a+1);i++) temp*=weightbwi(i);
      temp*=1/sqrt(pow(M_PI, p1a+1));
      bwi2 = xsmatrix.row(db);
      temp*=exp(-pow(rii.norm(), 2)/2)*exp(pow(bwi2.norm(), 2));
      temp*=getdeterminant(Hi2);
      temp/=sqrt(getdeterminant(Sig));
      dem+=temp;
    }
    loglike+=log(dem);
  }
  return loglike;
}