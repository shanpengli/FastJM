#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List updateH0(const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                    const Eigen::MatrixXd & X2, const Eigen::VectorXd & survtime, 
                    const Eigen::VectorXd & cmprsk,
                    const Eigen::MatrixXd & FUNE,
                    const Eigen::VectorXd & vi,
                    Eigen::MatrixXd & H01, Eigen::MatrixXd & H02){
  
  int a = H01.rows();
  int b = H02.rows();
  int k = X2.rows();
  
  int i,p,q,j,t,u;
  
  double scalef=0;
  
  /* calculate H01, H02*/
  double dem=0;
  int risk1_index=a-1;
  int risk2_index=b-1;
  
  for (j=0;j<k;j++)
  {
    
    dem+=FUNE(0,j)*exp(MultVV(X2.row(j), gamma1))*vi(j);
    
    if (cmprsk(j) == 1)
    {
      
      if (j == k-1)
      {
        H01(risk1_index, 2)=H01(risk1_index, 1)/dem;
        risk1_index--;
      }
      else if (survtime(j+1) != survtime(j))
      {
        H01(risk1_index, 2)=H01(risk1_index, 1)/dem;
        risk1_index--;
      }
      
      else
      {
        for (j=j+1;j<k;j++)
        {
          
          dem+=FUNE(0,j)*exp(MultVV(X2.row(j), gamma1))*vi(j);
          
          if (j == k-1)
          {
            H01(risk1_index, 2)=H01(risk1_index, 1)/dem;
            risk1_index--;
            break;
          }
          else if (survtime(j+1) != survtime(j))
          {
            H01(risk1_index, 2)=H01(risk1_index, 1)/dem;
            risk1_index--;
            break;
          }
          else continue;
        }
      }
      
    }
    else continue;
  }
  dem=0;
  for (j=0;j<k;j++)
  {
    
    dem+=FUNE(1,j)*exp(MultVV(X2.row(j), gamma2))*vi(j);
    
    if (cmprsk(j) == 2)
    {
      
      if (j == k-1)
      {
        H02(risk2_index, 2)=H02(risk2_index, 1)/dem;
        risk2_index--;
      }
      else if (survtime(j+1) != survtime(j))
      {
        H02(risk2_index, 2)=H02(risk2_index, 1)/dem;
        risk2_index--;
      }
      
      else
      {
        for (j=j+1;j<k;j++)
        {
          
          dem+=FUNE(1,j)*exp(MultVV(X2.row(j), gamma2))*vi(j);
          
          if (j == k-1)
          {
            H02(risk2_index, 2)=H02(risk2_index, 1)/dem;
            risk2_index--;
            break;
          }
          else if (survtime(j+1) != survtime(j))
          {
            H02(risk2_index, 2)=H02(risk2_index, 1)/dem;
            risk2_index--;
            break;
          }
          else continue;
        }
      }
      
    }
    else continue;
  }
  
  return Rcpp::List::create(Rcpp::Named("H01")= H01,
                            Rcpp::Named("H02")= H02);
}