// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include "basics.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

// another simple example: outer product of a vector,
// returning a matrix
//
// [[Rcpp::export]]
double  MultVV(const Eigen::VectorXd & x, const Eigen::VectorXd & y) {
    double v = x.transpose() * y;
    return v;
}

// [[Rcpp::export]]
Eigen::VectorXd CumSum(const Eigen::VectorXd & x) {
    Eigen::VectorXd a(x.size());
    int i;
    double temp = 0;
    for (i=0;i<a.size();i++) {
        temp += x(i);
        a(i) = temp;
    }
    return a;
}

// [[Rcpp::export]]
Eigen::MatrixXd MultVVoutprod(const Eigen::VectorXd & x) {
    Eigen::MatrixXd m = x * x.transpose();
    return m;
}

// [[Rcpp::export]]
Eigen::MatrixXd MultVV2outprod(const Eigen::VectorXd & x, const Eigen::VectorXd & y) {
    Eigen::MatrixXd m = x * y.transpose();
    return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double  MultVVinprod(const Eigen::VectorXd & x) {
    double v = x.transpose() * x;
    return v;
}

// [[Rcpp::export]]
double getdeterminant(const Eigen::MatrixXd & H) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd eigenSQ = svd.singularValues();
    double v=1;
    int i;
    for (i=0;i<eigenSQ.size();i++) v*=eigenSQ(i);
    return v;
}


// [[Rcpp::export]]
double  CH(const Eigen::MatrixXd & H, double t) {
    
    int a = H.rows();
    int i;
    double ch;
    if (t < H(0,0)) ch=0;
    else {
        ch=0;
        i=0;
        do {
            ch+=H(i, 2);
            i+=1;
            
        } while (i<=a-1 && t>= H(i,0));
        }
    
    return (ch);
    
}

// [[Rcpp::export]]
double HAZ(const Eigen::MatrixXd & H, double t) {
    
    int a = H.rows();
    int i;
    double temp=0;
    for (i=0;i<a;i++) {
        if (t == H(i, 0)) temp = H(i,2);
        }
    
    return(temp);
    
    }

// [[Rcpp::export]]
Eigen::MatrixXd MultMM(const Eigen::MatrixXd & x, const Eigen::MatrixXd & y) {
    Eigen::MatrixXd v = x * y;
    return v;
}

// [[Rcpp::export]]
double GetCIF1CR(const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                 const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                 const Eigen::VectorXd & X2,
                 const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                 const double s, const double u, const Eigen::VectorXd & bi){
  
  int a = H01.rows();
  int b = H02.rows();
  
  int i = 0;
  Eigen::VectorXd CH01 = Eigen::VectorXd::Zero(a);
  
  double hazard=0;
  for (i=0;i<a;i++) {
    hazard = hazard + H01(i, 2);
    CH01(i) = hazard;
  }
  
  Eigen::VectorXd FCH02 = Eigen::VectorXd::Zero(a);
  int count = 0;
  i = 0;
  while (count < b && i < a) {
    if (H02(count, 0) < H01(i, 0)) {
      FCH02(i) = FCH02(i) + H02(count, 2);
      count++;
    } else {
      i++;
    }
  }
  Eigen::VectorXd CH02 = Eigen::VectorXd::Zero(a);
  hazard=0;
  for (i=0;i<a;i++) {
    hazard = hazard + FCH02(i);
    CH02(i) = hazard;
  }
  
  double CIF1=0;
  
  for (i=0;i<a;i++) {
    if (s < H01(i, 0) && u >= H01(i, 0)) {
      if (i >= 1) {
        CIF1 = CIF1 + H01(i, 2)*exp(MultVV(X2,gamma1)+MultVV(alpha1,bi))*
          exp(-CH01(i-1)*exp(MultVV(X2,gamma1)+MultVV(alpha1,bi))-
          CH02(i-1)*exp(MultVV(X2,gamma2)+MultVV(alpha2,bi)));
      } else {
        CIF1 = CIF1 + H01(i, 2)*exp(MultVV(X2,gamma1)+MultVV(alpha1,bi));
      }
    } else continue;
    
  }
  
  return CIF1;
}

// [[Rcpp::export]]
double GetCIF2CR(const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                 const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                 const Eigen::VectorXd & X2,
                 const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                 const double s, const double u, const Eigen::VectorXd & bi){
  
  int a = H01.rows();
  int b = H02.rows();
  
  int i = 0;
  Eigen::VectorXd CH02 = Eigen::VectorXd::Zero(b);
  
  double hazard=0;
  for (i=0;i<b;i++) {
    hazard = hazard + H02(i, 2);
    CH02(i) = hazard;
  }
  
  Eigen::VectorXd FCH01 = Eigen::VectorXd::Zero(b);
  int count = 0;
  i = 0;
  while (count < a && i < b) {
    if (H01(count, 0) < H02(i, 0)) {
      FCH01(i) = FCH01(i) + H01(count, 2);
      count++;
    } else {
      i++;
    }
  }
  Eigen::VectorXd CH01 = Eigen::VectorXd::Zero(b);
  hazard=0;
  for (i=0;i<b;i++) {
    hazard = hazard + FCH01(i);
    CH01(i) = hazard;
  }
  
  double CIF2=0;
  
  for (i=0;i<b;i++) {
    if (s < H02(i, 0) && u >= H02(i, 0)) {
      if (i >= 1) {
        CIF2 = CIF2 + H02(i, 2)*exp(MultVV(X2,gamma2)+MultVV(alpha2,bi))*
          exp(-CH01(i-1)*exp(MultVV(X2,gamma1)+MultVV(alpha1,bi))-
          CH02(i-1)*exp(MultVV(X2,gamma2)+MultVV(alpha2,bi)));
      } else {
        CIF2 = CIF2 + H02(i, 2)*exp(MultVV(X2,gamma2)+MultVV(alpha2,bi));
      }
    } else continue;
    
  }
  
  return CIF2;
}

// [[Rcpp::export]]
Eigen::VectorXd GetCIF1CRall(const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                             const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                             const Eigen::VectorXd & X2,
                             const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                             const double s, const Eigen::VectorXd & timecif, const Eigen::VectorXd & bi){
  
  int a = H01.rows();
  int b = H02.rows();
  
  int i = 0;
  Eigen::VectorXd CH01 = Eigen::VectorXd::Zero(a);
  
  double hazard=0;
  for (i=0;i<a;i++) {
    hazard = hazard + H01(i, 2);
    CH01(i) = hazard;
  }
  
  Eigen::VectorXd FCH02 = Eigen::VectorXd::Zero(a);
  int count = 0;
  i = 0;
  while (count < b && i < a) {
    if (H02(count, 0) < H01(i, 0)) {
      FCH02(i) = FCH02(i) + H02(count, 2);
      count++;
    } else {
      i++;
    }
  }
  Eigen::VectorXd CH02 = Eigen::VectorXd::Zero(a);
  hazard=0;
  for (i=0;i<a;i++) {
    hazard = hazard + FCH02(i);
    CH02(i) = hazard;
  }
  
  Eigen::VectorXd CIF1 = Eigen::VectorXd::Zero(timecif.size());
  double CIF1temp=0;
  int countCIF1=0;
  
  for (i=0;i<a;i++) {
    if (s < H01(i, 0)) {
      if (i >= 1) {
        CIF1temp = CIF1temp + H01(i, 2)*exp(MultVV(X2,gamma1)+MultVV(alpha1,bi))*
          exp(-CH01(i-1)*exp(MultVV(X2,gamma1)+MultVV(alpha1,bi))-
          CH02(i-1)*exp(MultVV(X2,gamma2)+MultVV(alpha2,bi)));
      } else {
        CIF1temp = CIF1temp + H01(i, 2)*exp(MultVV(X2,gamma1)+MultVV(alpha1,bi));
      }
      CIF1(countCIF1)=CIF1temp;
      countCIF1++;
    } else continue;
  }
  
  return CIF1;
}

// [[Rcpp::export]]
Eigen::VectorXd GetCIF2CRall(const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                             const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                             const Eigen::VectorXd & X2,
                             const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                             const double s, const Eigen::VectorXd & timecif, const Eigen::VectorXd & bi){
  
  int a = H01.rows();
  int b = H02.rows();
  
  int i = 0;
  Eigen::VectorXd CH02 = Eigen::VectorXd::Zero(b);
  
  double hazard=0;
  for (i=0;i<b;i++) {
    hazard = hazard + H02(i, 2);
    CH02(i) = hazard;
  }
  
  Eigen::VectorXd FCH01 = Eigen::VectorXd::Zero(b);
  int count = 0;
  i = 0;
  while (count < a && i < b) {
    if (H01(count, 0) < H02(i, 0)) {
      FCH01(i) = FCH01(i) + H01(count, 2);
      count++;
    } else {
      i++;
    }
  }
  Eigen::VectorXd CH01 = Eigen::VectorXd::Zero(b);
  hazard=0;
  for (i=0;i<b;i++) {
    hazard = hazard + FCH01(i);
    CH01(i) = hazard;
  }
  
  Eigen::VectorXd CIF2 = Eigen::VectorXd::Zero(timecif.size());
  double CIF2temp=0;
  int countCIF2=0;
  
  for (i=0;i<b;i++) {
    if (s < H02(i, 0)) {
      if (i >= 1) {
        CIF2temp = CIF2temp + H02(i, 2)*exp(MultVV(X2,gamma2)+MultVV(alpha2,bi))*
          exp(-CH01(i-1)*exp(MultVV(X2,gamma1)+MultVV(alpha1,bi))-
          CH02(i-1)*exp(MultVV(X2,gamma2)+MultVV(alpha2,bi)));
      } else {
        CIF2temp = CIF2temp + H02(i, 2)*exp(MultVV(X2,gamma2)+MultVV(alpha2,bi));
      }
      CIF2(countCIF2)=CIF2temp;
      countCIF2++;
    } else continue;
    
  }
  
  return CIF2;
}


