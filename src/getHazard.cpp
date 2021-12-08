#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
int getHazard(const Eigen::Map<Eigen::VectorXd> & CumuH01,
               const Eigen::Map<Eigen::VectorXd> & CumuH02,
               const Eigen::Map<Eigen::VectorXd> & survtime,
               const Eigen::Map<Eigen::VectorXd> & cmprsk,
               const Eigen::Map<Eigen::MatrixXd> & H01,
               const Eigen::Map<Eigen::MatrixXd> & H02,
               Eigen::Map<Eigen::VectorXd> & CUH01,
               Eigen::Map<Eigen::VectorXd> & CUH02,
               Eigen::Map<Eigen::VectorXd> & HAZ01,
               Eigen::Map<Eigen::VectorXd> & HAZ02) {
  
  int a=H01.rows();
  int b=H02.rows();
  int k=survtime.size();
  
  int risk1_index=a-1;
  int risk2_index=b-1;
  
  int j;
  
  for (j=0;j<k;j++)
  {
    if (risk1_index>=0)
    {
      if (survtime(j) >= H01(risk1_index, 0))
      {
        CUH01(j) = CumuH01(risk1_index);

      }
      else
      {
        risk1_index--;
        if (risk1_index>=0)
        {
          CUH01(j) = CumuH01(risk1_index);
        }
      }
    }
    else
    {
      risk1_index=0;
    }
  }
  for (j=0;j<k;j++)
  {
    if (risk2_index>=0)
    {
      if (survtime(j) >= H02(risk2_index, 0))
      {
        CUH02(j) = CumuH02(risk2_index);
      }
      else
      {
        risk2_index--;
        if (risk2_index>=0)
        {
          CUH02(j) = CumuH02(risk2_index);
        }
      }
    }
    else
    {
      risk2_index=0;
    }
  }
  
  risk1_index = a-1;
  risk2_index = b-1;
  
  for (j=0;j<k;j++)
  {
    if (risk1_index>=0)
    {
      //gsl_vector_set(CUH01,j,gsl_vector_get(CumuH01,risk1_index));
      if (survtime(j) == H01(risk1_index, 0))
      {
        HAZ01(j) = H01(risk1_index, 2);
      }
      if (cmprsk(j) == 1)
      {
        if (j == k-1)
        {
          risk1_index--;
        }
        else if (survtime(j+1) != survtime(j))
        {
          risk1_index--;
        }
        else
        {
          for (j=j+1;j<k;j++)
          {
            if (survtime(j) == H01(risk1_index, 0))
            {
              HAZ01(j) = H01(risk1_index, 2);
            }
            if (j == k-1)
            {
              risk1_index--;
              break;
            }
            else if (survtime(j+1) != survtime(j))
            {
              risk1_index--;
              break;
            }
            else continue;
          }
        }
      }
      else continue;
    }
    else continue;
  }
  
  for (j=0;j<k;j++)
  {
    if (risk2_index>=0)
    {
      if (survtime(j) == H02(risk2_index, 0))
      {
        HAZ02(j) = H02(risk2_index, 2);
      }
      if (cmprsk(j) == 2)
      {
        if (j == k-1)
        {
          risk2_index--;
        }
        else if (survtime(j+1) != survtime(j))
        {
          risk2_index--;
        }
        else
        {
          for (j=j+1;j<k;j++)
          {
            if (survtime(j) == H02(risk2_index, 0))
            {
              HAZ02(j) = H02(risk2_index, 2);
            }
            
            if (j == k-1)
            {
              risk2_index--;
              break;
            }
            else if (survtime(j+1) != survtime(j))
            {
              risk2_index--;
              break;
            }
            else continue;
          }
        }
      }
      else continue;
    }
    else continue;
  }
  
  
  return 0;
  }