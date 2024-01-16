#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Eigen::MatrixXd getPusCR(const Eigen::VectorXd & beta, 
                         const Eigen::VectorXd & gamma1,
                         const Eigen::VectorXd & gamma2,  
                         const Eigen::VectorXd & alpha1,
                         const Eigen::VectorXd & alpha2,  
                         const Eigen::MatrixXd & Sig, const double sigma, 
                         const Eigen::MatrixXd & H01,
                         const Eigen::MatrixXd & H02,
                         const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
                         const Eigen::VectorXd & Y, const Eigen::MatrixXd & X2,
                         const Eigen::VectorXd & mdata, const Eigen::VectorXi & mdataS,
                         const Eigen::VectorXd & survtime,
                         const Eigen::MatrixXd & xsmatrix, 
                         const Eigen::MatrixXd & wsmatrix,
                         const double u) {
  
  int a = H01.rows();
  int b = H02.rows();
  //Calculate cumulative incidence function for type 1 failure
  Eigen::VectorXd CH011 = CumSum(H01.col(2));
  Eigen::VectorXd CH012 = Eigen::VectorXd::Zero(a);
  int count = 0;
  int i = 0;
  while ((count < b) && (i < a)) {
    if (H02(count, 0) <= H01(i, 0)) {
      CH012(i) += H02(count, 2);
      count += 1;
    } else {
      i += 1;
    }
  }
  CH012 = CumSum(CH012);
  
  Eigen::VectorXd CH021 = Eigen::VectorXd::Zero(b);
  Eigen::VectorXd CH022 = CumSum(H02.col(2));
  count = 0;
  i = 0;
  while ((count < a) && (i < b)) {
    if (H01(count, 0) <= H02(i, 0)) {
      CH021(i) += H01(count, 2);
      count += 1;
    } else {
      i += 1;
    }
  }
  CH021 = CumSum(CH021);
  
  double CIF1 = 0;
  double CIF2 = 0;
  
  //calculate the square root of random effect covariance matrix 
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd eigenSQ = svd.singularValues();
  int j,q,db;
  for (i=0;i<eigenSQ.size();i++) {
    eigenSQ(i) = sqrt(eigenSQ(i));
  }
  Eigen::MatrixXd SigSQRT  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  int p1a=Z.cols();
  int k=X2.rows();
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd weightbi(p1a);
  Eigen::VectorXd ri(p1a);
  Eigen::MatrixXd Pus(k, 2);
  
  int point=wsmatrix.rows();
  
  double dems,demu1,demu2,cuh01,cuh02,xgamma1,xgamma2,xalpha1,xalpha2,tempu1,tempu2,temps,mu,zb,temp;
  
  for (j=0;j<k;j++) {
    dems=0;
    demu1=0;
    demu2=0;
    q=mdata(j);
    cuh01=CH(H01, survtime(j));
    cuh02=CH(H02, survtime(j));
    xgamma1=MultVV(X2.row(j),gamma1);
    xgamma2=MultVV(X2.row(j),gamma2);
    
    for (db=0;db<point;db++) {
      CIF1=0;
      CIF2=0;
      bi = xsmatrix.row(db);
      weightbi = wsmatrix.row(db);
      ri = sqrt(2)*SigSQRT*bi;
      temp=exp(10);
      
      for (i=0;i<p1a;i++) bi(i)=ri(i);
      
      xalpha1=MultVV(bi, alpha1);
      xalpha2=MultVV(bi, alpha2);
      
      for (i=0;i<q;i++) {
        mu=MultVV(X1.row(mdataS(j)-1+i),beta);
        zb=MultVV(Z.row(mdataS(j)-1+i),bi);
        temp*=exp(-1/(2*sigma)*pow((Y(mdataS(j)-1+i) - mu - zb), 2)); 
      }
      
      for (i=0;i<p1a;i++) temp*=weightbi(i);
      
      ////Calculate cumulative incidence function for type 1 failure
      for (i=0;i<a;i++) {
        if ((survtime(j) <= H01(i, 0)) && (u >= H01(i, 0))) {
          CIF1 += H01(i, 2)*exp(xgamma1 + xalpha1)*
            exp(-CH011(i)*exp(xgamma1 + xalpha1)-
            CH012(i)*exp(xgamma2 + xalpha2));
        } else continue;
      }
      
      ////Calculate cumulative incidence function for type 2 failure
      for (i=0;i<b;i++) {
        if ((survtime(j) <= H02(i, 0)) && (u >= H02(i, 0))) {
          CIF2 += H02(i, 2)*exp(xgamma2 + xalpha2)*
            exp(-CH021(i)*exp(xgamma1 + xalpha1)-
            CH022(i)*exp(xgamma2 + xalpha2));
        } else continue;
      }
      
      temps=temp*exp(0-cuh01*exp(xgamma1+xalpha1)-cuh02*exp(xgamma2+xalpha2));
      tempu1=temp*CIF1;
      tempu2=temp*CIF2;
      
      dems+=temps;
      demu1+=tempu1;
      demu2+=tempu2;
    }
    Pus(j,0) = demu1/dems;
    Pus(j,1) = demu2/dems;
    
  }
  return(Pus);
  
  
  
}
