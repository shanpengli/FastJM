#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List getEC(const Eigen::VectorXd & beta, const Eigen::VectorXd & tau, 
           const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
           const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
           const double vee1, const double vee2, const Eigen::MatrixXd & H01, 
           const Eigen::MatrixXd & H02, const Eigen::MatrixXd & Sig, 
           const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
           const Eigen::MatrixXd & W, const Eigen::VectorXd & Y, 
           const Eigen::MatrixXd & X2, const Eigen::VectorXd & survtime, 
           const Eigen::VectorXd & cmprsk, const Eigen::VectorXd & mdata, 
           const Eigen::VectorXi & mdataS, const Eigen::MatrixXd & xsmatrix, 
           const Eigen::MatrixXd & wsmatrix, 
           const Eigen::VectorXd & CUH01, 
           const Eigen::VectorXd & CUH02, 
           const Eigen::VectorXd & HAZ01, 
           const Eigen::VectorXd & HAZ02){ 
  
  //calculate the square root of random effect covariance matrix 
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd eigenSQ = svd.singularValues();
  int i,j,q,t,db,u;
  for (i=0;i<eigenSQ.size();i++) {
    eigenSQ(i) = sqrt(eigenSQ(i));
    }
  Eigen::MatrixXd SigSQRT  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  int k=mdata.size();
  int p1a=Z.cols();
  
  double dem,cuh01,cuh02,haz01,haz02,xgamma1,xgamma2,temp,mu,sigma,zb,wi;
  Eigen::VectorXd bwi(p1a+1);
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd weightbwi(p1a+1);
  Eigen::VectorXd ri(p1a+1);
  
  /* Define functions of bi wi*/
  //exp(-w)
  Eigen::VectorXd FUNENW = Eigen::VectorXd::Zero(k);
  //bexp(-w)
  Eigen::MatrixXd FUNBENW = Eigen::MatrixXd::Zero(p1a,k);
  //b
  Eigen::MatrixXd FUNB = Eigen::MatrixXd::Zero(p1a,k);
  //w
  Eigen::VectorXd FUNW = Eigen::VectorXd::Zero(k);
  //bbT
  Eigen::MatrixXd FUNBS = Eigen::MatrixXd::Zero(p1a*(p1a+1)/2,k);
  //bw
  Eigen::MatrixXd FUNBW = Eigen::MatrixXd::Zero(p1a,k);
  //w2
  Eigen::VectorXd FUNWS = Eigen::VectorXd::Zero(k);
  //bbTexp(-w)
  Eigen::MatrixXd FUNBSENW = Eigen::MatrixXd::Zero(p1a*(p1a+1)/2,k);
  //exp(alpha*b + vee*w)
  Eigen::MatrixXd FUNEC = Eigen::MatrixXd::Zero(2,k);
  //bexp(alpha*b + vee*w)
  Eigen::MatrixXd FUNBEC = Eigen::MatrixXd::Zero(2*p1a,k);
  //bbTexp(alpha*b + vee*w)
  Eigen::MatrixXd FUNBSEC = Eigen::MatrixXd::Zero(2*p1a*(p1a+1)/2,k);
  //wexp(alpha*b + vee*w)
  Eigen::MatrixXd FUNWEC = Eigen::MatrixXd::Zero(2,k);
  //w2exp(alpha*b + vee*w)
  Eigen::MatrixXd FUNWSEC = Eigen::MatrixXd::Zero(2,k);
  
  int point=wsmatrix.rows();
  
  for(j=0;j<k;j++)
  {
    dem=0;
    q=mdata(j);
    cuh01=CUH01(j);
    cuh02=CUH02(j);
    haz01=HAZ01(j);
    haz02=HAZ02(j);
    
    // Rprintf("%d th subject cuh01 is: %f\n", j, cuh01);
    // 
    // cuh01=CH(H01, survtime(j));
    // cuh02=CH(H02, survtime(j));
    // haz01=HAZ(H01, survtime(j));
    // haz02=HAZ(H02, survtime(j));
    // 
    // Rprintf("%d th subject cuh01 is: %f\n", j, cuh01);
    // 
    xgamma1=MultVV(X2.row(j),gamma1);
    
    xgamma2=MultVV(X2.row(j),gamma2);
    
    for (db=0;db<point;db++) {
      
      bwi = xsmatrix.row(db);
      weightbwi = wsmatrix.row(db);
      ri = sqrt(2)*SigSQRT*bwi;
      temp=exp(10);
      
      // if (db<1) {
      //   Rprintf("%d th subject temp is: %f\n", j, temp);
      // }
      
      for (i=0;i<p1a;i++) bi(i)=ri(i);
      wi=ri(p1a);
      
      for (i=0;i<q;i++) {
        mu=MultVV(X1.row(mdataS(j)-1+i),beta);
        zb=MultVV(Z.row(mdataS(j)-1+i),bi);
        sigma=exp(MultVV(W.row(mdataS(j)-1+i),tau) + wi);
        temp*=1/sqrt(sigma)*exp(-1/(2*sigma)*pow((Y(mdataS(j)-1+i) - mu - zb), 2)); 
        }
      
      if(cmprsk(j)==1)  temp*=haz01*exp(xgamma1+MultVV(alpha1,bi)+vee1*wi);
      if(cmprsk(j)==2)  temp*=haz02*exp(xgamma2+MultVV(alpha2,bi)+vee2*wi);
      
      temp*=exp(0-cuh01*exp(xgamma1+MultVV(alpha1,bi)+vee1*wi)-cuh02*exp(xgamma2+MultVV(alpha2,bi)+vee2*wi));
      for (i=0;i<(p1a+1);i++) temp*=weightbwi(i);
      
      dem+=temp;
      
      //calculate h(bi,wi)
      FUNB.col(j)+=temp*bi;
      
      FUNW(j)+=wi*temp;
      FUNENW(j)+=exp(-wi)*temp;
      FUNBENW.col(j)+=exp(-wi)*temp*bi;
      for (i=0;i<p1a;i++) {
        FUNBS(i,j)+=temp*pow(bi(i),2);
        FUNBSENW(i,j)+=temp*exp(-wi)*pow(bi(i),2);
        }
      
      if (p1a > 1) {
        u=0;
        for(i=1;i<p1a;i++)
        {
          for(t=0;t<p1a-i;t++) {
            FUNBS(p1a+u,j) += temp*bi(t)*bi(t+i);
            FUNBSENW(p1a+u,j) += temp*exp(-wi)*bi(t)*bi(t+i);
            u++;
          }   
        }
      }
      
      FUNBW.col(j)+=wi*temp*bi;
      FUNWS(j)+=pow(wi,2)*temp;
      
      FUNEC(0,j)+=temp*exp(MultVV(alpha1,bi)+vee1*wi);
      FUNEC(1,j)+=temp*exp(MultVV(alpha2,bi)+vee2*wi);
      
      for (i=0;i<p1a;i++) {
        FUNBEC(i,j)+=temp*bi(i)*exp(MultVV(alpha1,bi)+vee1*wi);
        FUNBEC(p1a+i,j)+=temp*bi(i)*exp(MultVV(alpha2,bi)+vee2*wi);
        }
      
      for (i=0;i<p1a;i++) {
        FUNBSEC(i,j)+=temp*exp(MultVV(alpha1,bi)+vee1*wi)*pow(bi(i),2);
        FUNBSEC(p1a*(p1a+1)/2+i,j)+=temp*exp(MultVV(alpha2,bi)+vee2*wi)*pow(bi(i),2);
        }
      
      if (p1a > 1) {
        u=0;
        for(i=1;i<p1a;i++)
        {
          for(t=0;t<p1a-i;t++)
          {
            FUNBSEC(p1a+u,j)+=temp*exp(MultVV(alpha1,bi)+vee1*wi)*bi(t)*bi(t+i);
            FUNBSEC(p1a*(p1a+1)/2+p1a+u,j)+=temp*exp(MultVV(alpha2,bi)+vee2*wi)*bi(t)*bi(t+i);
            u++;
          }
        }
      }
      

      
      FUNWEC(0,j)+=temp*wi*exp(MultVV(alpha1,bi)+vee1*wi);
      FUNWEC(1,j)+=temp*wi*exp(MultVV(alpha2,bi)+vee2*wi);
      
      FUNWSEC(0,j)+=temp*pow(wi,2)*exp(MultVV(alpha1,bi)+vee1*wi);
      FUNWSEC(1,j)+=temp*pow(wi,2)*exp(MultVV(alpha2,bi)+vee2*wi);
      
      }
      
    if(dem==0) {
      Rprintf("E step ran into issue for the %dth subject. Program stops.\n", j);
      return ( 100.0 );
    } 
    
    FUNB.col(j)/=dem;
    FUNW(j)/=dem;
    FUNENW(j)/=dem;
    FUNBENW.col(j)/=dem;
    FUNBW.col(j)/=dem;
    FUNWS(j)/=dem;
    FUNBSENW.col(j)/=dem;
    FUNBS.col(j)/=dem;
    FUNEC.col(j)/=dem;
    FUNWEC.col(j)/=dem;
    FUNWSEC.col(j)/=dem;
    FUNBEC.col(j)/=dem;
    FUNBSEC.col(j)/=dem;
    // for (i=0;i<p1a;i++) {
    //   //FUNBENW(i,j)/=dem;
    //   //FUNBW(i,j)/=dem;
    //   }
    // 
    // for (i=0;i<(p1a*(p1a+1)/2);i++) {
    //   //FUNBSENW(i,j)/=dem;
    //   //FUNBS(i,j)/=dem;
    //   }
    // 
    // 
    // 
    // for (i=0;i<2;i++) {
    //   FUNEC(i,j)/=dem;
    //   FUNWEC(i,j)/=dem;
    //   FUNWSEC(i,j)/=dem;
    //   }
    // 
    // for(i=0;i<(2*p1a);i++) {
    //   FUNBEC(i,j)/=dem;
    //   }
    // 
    // for (i=0;i<(2*p1a*(p1a+1)/2);i++) {
    //   FUNBSEC(i,j)/=dem;
    //   }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("FUNB")=FUNB,
                            Rcpp::Named("FUNW")=FUNW,
                            Rcpp::Named("FUNENW")=FUNENW,
                            Rcpp::Named("FUNBENW")=FUNBENW,
                            Rcpp::Named("FUNBW")=FUNBW,
                            Rcpp::Named("FUNWS")=FUNWS,
                            Rcpp::Named("FUNBSENW")=FUNBSENW,
                            Rcpp::Named("FUNBS")=FUNBS,
                            Rcpp::Named("FUNEC")=FUNEC,
                            Rcpp::Named("FUNBEC")=FUNBEC,
                            Rcpp::Named("FUNBSEC")=FUNBSEC,
                            Rcpp::Named("FUNWEC")=FUNWEC,
                            Rcpp::Named("FUNWSEC")=FUNWSEC);
}


