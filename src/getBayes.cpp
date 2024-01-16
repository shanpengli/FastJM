#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List getBayes(const Eigen::VectorXd & beta,
                  const Eigen::MatrixXd & Sig, 
                  const double sigma, const Eigen::MatrixXd & Z, 
                  const Eigen::MatrixXd & X1, const Eigen::VectorXd & Y, 
                  const Eigen::VectorXd & mdata, 
                  const Eigen::VectorXi & mdataS) {
  
  int p1a = Sig.cols();
  int k = mdata.size();
  int p1 = beta.size();
  
  Eigen::MatrixXd Posbi = Eigen::MatrixXd::Zero(k, p1a);
  Eigen::MatrixXd Poscov = Eigen::MatrixXd::Zero(k*p1a, p1a);
  Eigen::MatrixXd XVX = Eigen::MatrixXd::Zero(p1, p1);
  Eigen::MatrixXd Hi = Eigen::MatrixXd::Zero(p1a, p1a);
  
  int q,i,j;
  
  for (j=0;j<k;j++) {

    q = mdata(j);

    Eigen::MatrixXd X1tilde = Eigen::MatrixXd::Zero(q, p1a);
    Eigen::MatrixXd X11 = Eigen::MatrixXd::Zero(q, p1);
    Eigen::VectorXd Yi = Eigen::VectorXd::Zero(q);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(q, q);

    for (i=0;i<q;i++) {
      X1tilde.row(i) = Z.row(mdataS(j)-1+i);
      X11.row(i) = X1.row(mdataS(j)-1+i);
      Yi(i) = Y(mdataS(j)-1+i);
    }
    
    V = X1tilde*Sig*X1tilde.transpose();

    for (i=0;i<q;i++) V(i,i) += sigma;

    XVX += X11.transpose()*V.inverse()*X11;

    Posbi.row(j) = Sig*X1tilde.transpose()*V.inverse()*(Yi - X11*beta);
  }

  for (j=0;j<k;j++) {

    q = mdata(j);
    Eigen::MatrixXd X1tilde = Eigen::MatrixXd::Zero(q, p1a);
    Eigen::MatrixXd X11 = Eigen::MatrixXd::Zero(q, p1);
    Eigen::VectorXd Yi = Eigen::VectorXd::Zero(q);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(q, q);

    for (i=0;i<q;i++) {
      X1tilde.row(i) = Z.row(mdataS(j)-1+i);
      X11.row(i) = X1.row(mdataS(j)-1+i);
      Yi(i) = Y(mdataS(j)-1+i);
    }

    V = X1tilde*Sig*X1tilde.transpose();
    
    for (i=0;i<q;i++) V(i,i) += sigma;

    Hi = Sig - Sig*X1tilde.transpose()*(V.inverse() -
      V.inverse()*X11*XVX.inverse()*X11.transpose()* V.inverse())*X1tilde*Sig;

    for (i=0;i<p1a;i++) Poscov.row(j*p1a+i) = Hi.row(i);

  }
  
  return Rcpp::List::create(Rcpp::Named("Posbi")=Posbi,
                            Rcpp::Named("Poscov")=Poscov);
  
  
}
