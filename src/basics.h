#include <RcppEigen.h>
#include <stdio.h>
#include <chrono>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

// [[Rcpp::depends(RcppEigen)]]

double  MultVV(const Eigen::VectorXd & x, const Eigen::VectorXd & y);

Eigen::MatrixXd MultVVoutprod(const Eigen::VectorXd & x);

Eigen::VectorXd CumSum(const Eigen::VectorXd & x);

double MultVVinprod(const Eigen::VectorXd & x);

double CH(const Eigen::MatrixXd & H, double t);

double HAZ(const Eigen::MatrixXd & H, double t);

Eigen::MatrixXd MultVV2outprod(const Eigen::VectorXd & x, const Eigen::VectorXd & y);

double getdeterminant(const Eigen::MatrixXd & H);
