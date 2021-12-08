#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List GetrisksetCSF(const Eigen::MatrixXd & cdata) {
  
  int k = cdata.rows();
  int u=0,a=0;
  int i,j;
  Eigen::MatrixXd FH01 = Eigen::MatrixXd::Zero(k, 2);
  
  /* find # events for risk 1*/
  for (j=0;j<k;j++)
  {
    if (cdata(j,1) == 1)
    {
      u++;
      if (j == k-1)
      {
        a++;
        FH01(k-a,0) = cdata(j,0);
        FH01(k-a,1) = u;
        u=0;
      }
      else if (cdata(j+1,0) != cdata(j,0))
      {
        a++;
        FH01(k-a,0) = cdata(j,0);
        FH01(k-a,1) = u;
        u=0;
      }
      else
      {
        for (j=j+1;j<k;j++)
        {
          if (cdata(j,1) == 1)
          {
            u++;
            if (j == k-1)
            {
              a++;
              FH01(k-a,0) = cdata(j,0);
              FH01(k-a,1) = u;
              u=0;
              break;
            }
            else if (cdata(j+1,0) != cdata(j,0))
            {
              a++;
              FH01(k-a,0) = cdata(j,0);
              FH01(k-a,1) = u;
              u=0;
              break;
            }
            else continue;
          }
          else
          {
            if (j == k-1)
            {
              a++;
              FH01(k-a,0) = cdata(j,0);
              FH01(k-a,1) = u;
              u=0;
              break;
            }
            else if (cdata(j+1,0) != cdata(j,0))
            {
              a++;
              FH01(k-a,0) = cdata(j,0);
              FH01(k-a,1) = u;
              u=0;
              break;
            }
            else continue;
          }
        }
      }
      
    }
    else continue;
  }
  
  if(a==0)
  {
    Rprintf("No failure time information for risk 1; Program exits\n");
    return ( -1.0 );
  } 
  
  Eigen::MatrixXd H01 = Eigen::MatrixXd::Zero(a, 3);
  for(i=0;i<3;i++)
  {
    if(i<=1)
    {
      for(j=a;j>0;j--)    H01(a-j,i) = FH01(k-j,i);
    }
    if(i==2)
    {
      for(j=0;j<a;j++)    H01(j,i) = 0.0001;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("H01")=H01);
  
}