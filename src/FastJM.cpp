//
//  FastJM.cpp
//  FastJM
//
//  Created by Shanpeng Li on 6/16/20.
//

#include <stdio.h>
#include "jmcs.hpp"
#include "jmcsf.hpp"
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List  jmcs_main(SEXP k, SEXP n1,SEXP p1,SEXP p2, SEXP p1a, SEXP maxiter, SEXP point, SEXP xs,SEXP ws,SEXP yfile, SEXP cfile, SEXP mfile, SEXP Betasigmafile, SEXP Sigcovfile, SEXP trace)
  {
    Rcpp::List result;
  try {

     result=jmcsspace::jmcs_cmain( as<int> (k), as<int> (n1), as<int> (p1), as<int> (p2),
            as<int> (p1a), as<int>(maxiter),as<int>(point), as<std::vector<double> >(xs),as<std::vector<double> >(ws), as<std::string> (yfile),
           as<std::string> (cfile),as<std::string>(mfile),as<std::string>(Betasigmafile),as<std::string>(Sigcovfile),as<int> (trace));
   if(Rf_isNull(result)){
        throw std::range_error("Possible files reading or format errors");
        }
     return result;
    } catch(std::exception &ex) {
    forward_exception_to_r(ex);
    } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;             // not reached

}

// [[Rcpp::export]]
Rcpp::List  jmcsf_main(SEXP k, SEXP n1,SEXP p1,SEXP p2, SEXP p1a, SEXP maxiter, SEXP point, SEXP xs,SEXP ws,SEXP yfile, SEXP cfile, SEXP mfile, SEXP Betasigmafile, SEXP Sigcovfile, SEXP trace)
  {
    Rcpp::List result;
  try {

     result=jmcsfspace::jmcsf_cmain( as<int> (k), as<int> (n1), as<int> (p1), as<int> (p2),
            as<int> (p1a), as<int>(maxiter),as<int>(point), as<std::vector<double> >(xs),as<std::vector<double> >(ws), as<std::string> (yfile),
           as<std::string> (cfile),as<std::string>(mfile),as<std::string>(Betasigmafile),as<std::string>(Sigcovfile),as<int> (trace));
   if(Rf_isNull(result)){
        throw std::range_error("Possible files reading or format errors");
        }
     return result;
    } catch(std::exception &ex) {
    forward_exception_to_r(ex);
    } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;             // not reached

}
