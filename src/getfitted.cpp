#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List getfitted(const Eigen::VectorXd & beta, const Eigen::MatrixXd & Z, 
                     const Eigen::MatrixXd & X1, const Eigen::VectorXd & Y, 
                     const Eigen::VectorXd & mdata, 
                     const Eigen::VectorXi & mdataS,
                     const Eigen::MatrixXd & FUNB) {
  
  int k = mdata.size();
  
  Eigen::VectorXd fitted = Eigen::VectorXd::Zero(mdata.sum());
  Eigen::VectorXd fittedmar = Eigen::VectorXd::Zero(mdata.sum());
  Eigen::VectorXd resid = Eigen::VectorXd::Zero(mdata.sum());
  Eigen::VectorXd residmar = Eigen::VectorXd::Zero(mdata.sum());

  int q,i,j;
  
  for (j=0;j<k;j++) {
    q = mdata(j);
    for (i=0;i<q;i++) fitted(mdataS(j)-1+i) = MultVV(X1.row(mdataS(j)-1+i), beta) + MultVV(Z.row(mdataS(j)-1+i), FUNB.col(j));
  }
  
  fittedmar = X1*beta;
  resid = Y - fitted;
  residmar = Y - fittedmar;

  return Rcpp::List::create(Rcpp::Named("resid")=resid,
                            Rcpp::Named("fitted")=fitted,
                            Rcpp::Named("fittedmar")=fittedmar,
                            Rcpp::Named("residmar")=residmar);
  
}
