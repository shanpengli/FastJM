#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Eigen::MatrixXd getfittedSurvSF(const Eigen::VectorXd & gamma1, const Eigen::MatrixXd & X2,
                                const Eigen::MatrixXd & CH01, const Eigen::VectorXd & alpha1, 
                                const Eigen::MatrixXd & FUNB) {
  
  int k = X2.rows();
  int q = CH01.rows();
  Eigen::MatrixXd fittedSurv = Eigen::MatrixXd::Zero(q, 2);
  Eigen::VectorXd a = Eigen::VectorXd::Zero(q);
  fittedSurv.col(0) = CH01.col(0);
  int i, j;
  for (j=0;j<k;j++) {
    a = -exp(MultVV(X2.row(j), gamma1) + MultVV(alpha1, FUNB.col(j)))*CH01.col(3);
    for (i=0;i<q;i++) a(i) = exp(a(i));
    fittedSurv.col(1) += a;
    }
  fittedSurv.col(1)/=k;
  
  return (fittedSurv);
  
}
