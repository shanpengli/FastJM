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

double GetCIF1CR(const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                 const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                 const Eigen::VectorXd & X2,
                 const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                 const double s, const double u, const Eigen::VectorXd & bi);
  

double GetCIF2CR(const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                 const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                 const Eigen::VectorXd & X2,
                 const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                 const double s, const double u, const Eigen::VectorXd & bi);


  
