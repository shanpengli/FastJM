#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Eigen::MatrixXd getfittedSurv(const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                              const Eigen::MatrixXd & X2, Eigen::MatrixXd & CH012, 
                              const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                              const Eigen::MatrixXd & FUNB) {
  
  int k = X2.rows();
  int q = CH012.rows();
  Eigen::MatrixXd fittedSurv = Eigen::MatrixXd::Zero(q, 2);
  Eigen::VectorXd a = Eigen::VectorXd::Zero(q);
  fittedSurv.col(0) = CH012.col(0);
  int i,j;
  for (j=1;j<q;j++) {
    if (CH012(j, 1) == 0) CH012(j, 1) = CH012(j-1, 1);
    if (CH012(j, 2) == 0) CH012(j, 2) = CH012(j-1, 2);
  }
  
  for (j=0;j<k;j++) {
    a = -exp(MultVV(X2.row(j), gamma1) + MultVV(alpha1, FUNB.col(j)))*CH012.col(1)-
      exp(MultVV(X2.row(j), gamma2) + MultVV(alpha2, FUNB.col(j)))*CH012.col(2);
    for (i=0;i<q;i++) a(i) = exp(a(i));
    fittedSurv.col(1) += a;
  }
  fittedSurv.col(1)/=k;
  
  return (fittedSurv);
  
}
