#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List OLS(const Eigen::MatrixXd & X, const Eigen::VectorXd & Y) {
  
  Eigen::LLT<Eigen::MatrixXd> LLT_of_K(X.transpose() * X);
  if ( !LLT_of_K.info() ) {
    Eigen::VectorXd betahat = LLT_of_K.solve(X.transpose() * Y);
    Eigen::VectorXd resid = Y - X * betahat;
    return Rcpp::List::create(Rcpp::Named("betahat")=betahat,
                              Rcpp::Named("resid")=resid);
  } else {
    return ( -1.0 );
  }
}