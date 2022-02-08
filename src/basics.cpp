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


