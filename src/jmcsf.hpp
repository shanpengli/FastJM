//
//  jmcsf.hpp
//  FastJM
//
//  Created by Shanpeng Li on 6/16/20.
//

#ifndef jmcsf_hpp
#define jmcsf_hpp

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <Rcpp.h>
using namespace Rcpp;

namespace jmcsfspace {

    double HAZ(const gsl_matrix *H, const double t);

    double CH(const gsl_matrix *H, const double t);

    double MulVV(const gsl_vector *Z,const gsl_vector *beta);

    void MulV(const gsl_vector *Z,gsl_matrix *ZZ);


    void MulMM(const gsl_matrix *A,const gsl_matrix *B,gsl_matrix *AB);

    int inv_matrix(gsl_matrix *x_square);

    double Abs(const double a, const double b);

    int DiffM1(const gsl_matrix *matrixa, const gsl_matrix *matrixb);

    double DiffM(const gsl_matrix *matrixa, const gsl_matrix *matrixb);

    double DiffV(const gsl_vector *veca, const gsl_vector *vecb);

    double Min(const double t1, const double t2);

    void STAT(gsl_matrix *store,int i,double *mean,double *sd);

    int GetN(double t);

    int DiffM2(const gsl_matrix *preH1,const gsl_matrix *H1);

    void TransM(const gsl_matrix *A, gsl_matrix *B);

    int Sbeta(gsl_vector *beta, double *sigma, const gsl_matrix *Y, const int p1a);



    int EM(
           gsl_vector *beta,
           gsl_matrix *gamma,
           gsl_vector *vee1,
           gsl_matrix *H01,
           double *sigma,
           gsl_matrix *sig,
           const int p1a,
           const gsl_matrix *Y,
           const gsl_matrix *C,
           const gsl_vector *M1,
           const gsl_matrix *Posbi,
           const gsl_matrix *Poscov,
           const int point,
           const std::vector<double> xs,
           const std::vector<double> ws
           );


    int GetCov(
               gsl_matrix *Cov,
               const gsl_vector *beta,
               const gsl_matrix *gamma,
               const gsl_vector *vee1,
               const gsl_matrix *H01,
               const double sigma,
               const gsl_matrix *sig,
               const gsl_matrix *Y,
               const gsl_matrix *C,
               const gsl_vector *M1,
               const gsl_matrix *Posbi,
               const gsl_matrix *Poscov,
               const int p1a,
               const int point,
               const std::vector<double> xs,
               const std::vector<double> ws
               );


    int GetE(
             gsl_matrix *FUNB,
             gsl_matrix *FUNBS,
             gsl_matrix *FUNBSE,
             gsl_matrix *FUNBE,
             gsl_matrix *FUNE,
             const gsl_vector *beta,
             const gsl_matrix *gamma,
             const gsl_vector *vee1,
             const gsl_matrix *H01,
             const double sigma,
             const gsl_matrix *sig,
             const gsl_matrix *Y,
             const gsl_matrix *C,
             const gsl_vector *M1,
             const gsl_matrix *Posbi,
             const gsl_matrix *Poscov,
             const int p1a,
             const int point,
             const std::vector<double> xs,
             const std::vector<double> ws
             );

    double Getloglik(
                     const gsl_vector *beta,
                     const gsl_matrix *gamma,
                     const gsl_vector *vee1,
                     const gsl_matrix *H01,
                     const double sigma,
                     const gsl_matrix *sig,
                     const gsl_matrix *Y,
                     const gsl_matrix *C,
                     const gsl_vector *M1,
                     const gsl_matrix *Posbi,
                     const gsl_matrix *Poscov,
                     const int p1a,
                     const int point,
                     const std::vector<double> xs,
                     const std::vector<double> ws
                     );

    int Diff(
             const gsl_vector *prebeta,
             const gsl_vector *beta,
             const gsl_matrix *pregamma,
             const gsl_matrix *gamma,
             const gsl_vector *prevee1,
             const gsl_vector *vee1,
             const gsl_matrix *preH01,
             const gsl_matrix *H01,
             const double presigma,
             const double sigma,
             const gsl_matrix *presig,
             const gsl_matrix *sig
             );

    double GetPosbi(
                    const gsl_matrix *Y,
                    const gsl_vector *beta,
                    const gsl_vector *M1,
                    const gsl_matrix *Sigcov,
                    const double sigma,
                    gsl_matrix *Xinv,
                    gsl_matrix *Posbi,
                    const int k,
                    const int p1,
                    const int p1a
                    );

    double GetPoscov(
                     const gsl_matrix *Y,
                     const gsl_vector *beta,
                     const gsl_vector *M1,
                     const gsl_matrix *Sigcov,
                     const double sigma,
                     const gsl_matrix *Xinv,
                     gsl_matrix *Poscov,
                     const int k,
                     const int p1,
                     const int p1a
    );

    //declare before use it
    Rcpp::List jmcsf_cmain(int k, int n1,int p1,int p2, int p1a, int maxiter, int point,std::vector<double> xs,  std::vector<double> ws, std::string yfile, std::string cfile, std::string mfile, std::string Betasigmafile, std::string Sigcovfile, int trace);


}

#endif /* jmcsf_hpp */
