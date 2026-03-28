#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List getCovSF_JMH(const Eigen::VectorXd & beta, const Eigen::VectorXd & tau, 
                  const Eigen::VectorXd & gamma1, const Eigen::VectorXd & alpha1, 
                  const double vee1, const Eigen::MatrixXd & H01, const Eigen::MatrixXd & Sig, 
                  const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
                  const Eigen::MatrixXd & W, const Eigen::VectorXd & Y, 
                  const Eigen::MatrixXd & X2, const Eigen::VectorXd & survtime, 
                  const Eigen::VectorXd & cmprsk, const Eigen::VectorXd & mdata, 
                  const Eigen::VectorXi & mdataS,
                  const Eigen::VectorXd & FUNENW, 
                  const Eigen::MatrixXd & FUNBENW, 
                  const Eigen::MatrixXd & FUNBS, 
                  const Eigen::MatrixXd & FUNBW, 
                  const Eigen::VectorXd & FUNWS, 
                  const Eigen::MatrixXd & FUNBSENW, 
                  const Eigen::MatrixXd & FUNEC, 
                  const Eigen::MatrixXd & FUNBEC,
                  const Eigen::MatrixXd & FUNBSEC, 
                  const Eigen::MatrixXd & FUNWEC, 
                  const Eigen::MatrixXd & FUNWSEC,
                  const Eigen::MatrixXd & FUNB,
                  const Eigen::VectorXd & FUNW) {


  int d = beta.size() + W.cols() + gamma1.size() + alpha1.size() +
    1 + Sig.cols()*(Sig.cols() + 1)/2;
  int a = H01.rows();
  int k = mdata.size();
  
  int i,p,q,j,t,u;
  
  double temp,temp1;
  
  Eigen::MatrixXd SS  = Eigen::MatrixXd::Zero(d,d);
  Eigen::MatrixXd SSinv  = Eigen::MatrixXd::Zero(d,d);
  Eigen::VectorXd S = Eigen::VectorXd::Zero(d);
  
  int risk1_index;
  int risk1_index_temp=a-1;
  int risk1_index_ttemp=a-1;
  int risk1_index_tttemp=a-1;
  int risk1_index_vtemp=a-1;
  int risk1_index_vttemp=a-1;
  int risk1_index_vtttemp=a-1;
  int risk1_index_wtemp=a-1;
  int risk1_index_wttemp=a-1;
  int risk1_index_wtttemp=a-1;

  Eigen::VectorXd CumuH01 = Eigen::VectorXd::Zero(a);
  
  temp1=0;
  for (j=0;j<a;j++) {
    temp1 += H01(j, 2);
    CumuH01(j) = temp1;
  }
  
  double epsilon=0;
  double qq=0;
  int p1 = beta.size();
  int p1b = tau.size();
  int p1a = Z.cols();
  int p2 = gamma1.size();

  Eigen::VectorXd SZ = Eigen::VectorXd::Zero(p1);
  Eigen::VectorXd SZ1 = Eigen::VectorXd::Zero(p1);
  Eigen::MatrixXd SZZ = Eigen::MatrixXd::Zero(p1, p1a);
  Eigen::VectorXd SZtau = Eigen::VectorXd::Zero(p1b);
  Eigen::MatrixXd bbT  = Eigen::MatrixXd::Zero(p1a, p1a);
  Eigen::MatrixXd bbT2  = Eigen::MatrixXd::Zero(p1a, p1a);

  Eigen::VectorXd X = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SX = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SX1 = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SRX = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd SRXX = Eigen::VectorXd::Zero(p2);
  Eigen::MatrixXd SXX1 = Eigen::MatrixXd::Zero(p2, a);
  Eigen::MatrixXd SXX11 = Eigen::MatrixXd::Zero(p2, a);
  Eigen::MatrixXd SRXX1 = Eigen::MatrixXd::Zero(p2, k);

  Eigen::VectorXd N = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TN = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TN1 = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TRN = Eigen::VectorXd::Zero(p1a);
  Eigen::VectorXd TRNN = Eigen::VectorXd::Zero(p1a);
  Eigen::MatrixXd TRNN1 = Eigen::MatrixXd::Zero(p1a,k);
  Eigen::MatrixXd TNN1 = Eigen::MatrixXd::Zero(p1a,a);
  Eigen::MatrixXd TNN11 = Eigen::MatrixXd::Zero(p1a,a);

  double G = 0;
  double TG = 0;
  double TG1 = 0;
  double TRG = 0;
  double TRGG = 0;
  Eigen::VectorXd TRGG1 = Eigen::VectorXd::Zero(k);
  Eigen::VectorXd TGG1 = Eigen::VectorXd::Zero(a);
  Eigen::VectorXd TGG11 = Eigen::VectorXd::Zero(a);

  Eigen::MatrixXd bsw = Eigen::MatrixXd::Zero(p1a+1,p1a+1);
  Eigen::MatrixXd bsw2 = Eigen::MatrixXd::Zero(p1a+1,p1a+1);

  for (j=0;j<k;j++) {
  
    S = Eigen::VectorXd::Zero(d);
    q = mdata(j);
  
    /* calculate score for beta */
    SZ = Eigen::VectorXd::Zero(p1);
    SZ1 = Eigen::VectorXd::Zero(p1);
    for (i=0;i<q;i++) {
  
      SZ1 += (Y(mdataS(j)-1+i) - MultVV(X1.row(mdataS(j)-1+i), beta))*
        exp(-MultVV(W.row(mdataS(j)-1+i),tau))*FUNENW(j)*X1.row(mdataS(j)-1+i);
  
      SZ += exp(-MultVV(W.row(mdataS(j)-1+i),tau))*
        MultVV(Z.row(mdataS(j)-1+i), FUNBENW.col(j))*X1.row(mdataS(j)-1+i);
  
    }
    for (i=0;i<p1;i++) S(i) = SZ1(i) - SZ(i);
  
    /* calculate score for tau */
    SZtau = Eigen::VectorXd::Zero(p1b);
    q=mdata(j);

    for(t=0;t<p1a;t++)   bbT(t,t) = FUNBSENW(t,j);
    if(p1a>1)
    {
      u=0;
      for(i=1;i<p1a;i++)
      {
        for(t=0;t<p1a-i;t++) {
          bbT(t,i+t) = FUNBSENW(p1a+u,j);
          bbT(i+t,t) = bbT(t,i+t);
          u++;
        }
      }
    }

    for (i=0;i<q;i++) {

      epsilon = Y(mdataS(j)-1+i) - MultVV(X1.row(mdataS(j)-1+i), beta);
      bbT2 = MultVVoutprod(Z.row(mdataS(j)-1+i))*bbT;
      qq = pow(epsilon, 2)*FUNENW(j) - 2*epsilon*MultVV(Z.row(mdataS(j)-1+i), FUNBENW.col(j)) + bbT2.trace();
      SZtau += 0.5*(exp(-MultVV(W.row(mdataS(j)-1+i),tau))*qq-1)*W.row(mdataS(j)-1+i);
    }
    for (i=0;i<p1b;i++) S(p1+i) = SZtau(i);
    
    
    /* calculate score for gamma */
    if (j == 0)
    {
      temp=0;
      risk1_index=risk1_index_temp;
      for (q=j;q<k;q++)
      {
    
        temp+=exp(MultVV(X2.row(q), gamma1))*FUNEC(0,q);
        SX += exp(MultVV(X2.row(q), gamma1))*FUNEC(0,q)*X2.row(q);
    
        if (cmprsk(q) == 1)
        {
          if (q == k-1)
          {
            SX*= H01(risk1_index, 1)/pow(temp, 2);
            SXX1.col(a-1-risk1_index) = SX;
            SRXX += SX;
            SX/= H01(risk1_index, 1)/pow(temp, 2);
            SX/=temp;
            SXX11.col(a-1-risk1_index) = SX;
            SX*=temp;
            risk1_index--;
          }
          else if (survtime(q+1) != survtime(q))
          {
            SX*= H01(risk1_index, 1)/pow(temp, 2);
            SXX1.col(a-1-risk1_index) = SX;
            SRXX += SX;
            SX/= H01(risk1_index, 1)/pow(temp, 2);
            SX/=temp;
            SXX11.col(a-1-risk1_index) = SX;
            SX*=temp;
            risk1_index--;
          }
          else
          {
            for (q=q+1;q<k;q++)
            {
              temp+=exp(MultVV(X2.row(q), gamma1))*FUNEC(0,q);
              SX += exp(MultVV(X2.row(q), gamma1))*FUNEC(0,q)*X2.row(q);
    
              if (q == k-1)
              {
                SX*= H01(risk1_index, 1)/pow(temp, 2);
                SXX1.col(a-1-risk1_index) = SX;
                SRXX += SX;
                SX/= H01(risk1_index, 1)/pow(temp, 2);
                SX/=temp;
                SXX11.col(a-1-risk1_index) = SX;
                SX*=temp;
                risk1_index--;
                break;
              }
              else if (survtime(q+1) != survtime(q))
              {
                SX*= H01(risk1_index, 1)/pow(temp, 2);
                SXX1.col(a-1-risk1_index) = SX;
                SRXX += SX;
                SX/= H01(risk1_index, 1)/pow(temp, 2);
                SX/=temp;
                SXX11.col(a-1-risk1_index) = SX;
                SX*=temp;
                risk1_index--;
                break;
              }
              else continue;
            }
          }
    
        }
        else continue;
      }
      SRXX1.col(j) = SRXX;
    }
    else
    {
      if (risk1_index_temp>=0)
      {
        if (survtime(j) >= H01(risk1_index_temp, 0))
        {
          SRXX1.col(j) = SRXX1.col(j-1);
        }
        else
        {
          risk1_index_temp--;
          if (risk1_index_temp>=0)
          {
            SRXX = SRXX1.col(j-1);
            SRXX -= SXX1.col(a-1-risk1_index_temp-1);
            SRXX1.col(j) = SRXX;
          }
        }
      }
      else
      {
        risk1_index_temp=0;
      }
    }
    SRX = SRXX1.col(j);
    
    if (j==0)
    {
      SRX -= CumuH01(risk1_index_ttemp)*X2.row(j);
    }
    else if (survtime(j) >= H01(risk1_index_ttemp, 0))
    {
      SRX -= CumuH01(risk1_index_ttemp)*X2.row(j);
    }
    else
    {
      risk1_index_ttemp--;
      if (risk1_index_ttemp>=0)
      {
        SRX -= CumuH01(risk1_index_ttemp)*X2.row(j);
      }
      else
      {
        SRX = SRXX1.col(j);
        risk1_index_ttemp=0;
      }
    }
    
    SRX*= exp(MultVV(X2.row(j), gamma1))*FUNEC(0,j);
    
    if (survtime(j) >= H01(risk1_index_tttemp, 0))
    {
      if (cmprsk(j) == 1)
      {
        X = X2.row(j);
        X -= SXX11.col(a-1-risk1_index_tttemp);
        X += SRX;
        for (q=0;q<p2;q++) S(p1+p1b+q) = X(q);
      }
      else
      {
        for (q=0;q<p2;q++) S(p1+p1b+q) = SRX(q);
      }
    }
    else
    {
      risk1_index_tttemp--;
      if (risk1_index_tttemp>=0)
      {
        if (cmprsk(j) == 1)
        {
          X = X2.row(j);
          X -= SXX11.col(a-1-risk1_index_tttemp);
          X += SRX;
          for (q=0;q<p2;q++) S(p1+p1b+q) = X(q);
        }
        else
        {
          for (q=0;q<p2;q++) S(p1+p1b+q) = SRX(q);
        }
      }
      else
      {
        risk1_index_tttemp=0;
        for (q=0;q<p2;q++) S(p1+p1b+q) = SRX(q);
      }
    }

    
    /* calculate score for alpha */
    /*  alpha1 */
    if (j == 0)
    {
      temp=0;
    
      TN = Eigen::VectorXd::Zero(p1a);
      TRN = Eigen::VectorXd::Zero(p1a);
    
      risk1_index=risk1_index_vtemp;
      for (q=j;q<k;q++)
      {
        temp += exp(MultVV(X2.row(q), gamma1))*FUNEC(0, q);
        for (i=0;i<p1a;i++) N(i) = FUNBEC(i,q);
        TN += exp(MultVV(X2.row(q), gamma1))*N;
        if (cmprsk(q) == 1)
        {
          if (q == k-1)
          {
            TN *= H01(risk1_index, 1)/pow(temp,2);
            TNN1.col(a-1-risk1_index) = TN;
            TRNN += TN;
            TN /= H01(risk1_index, 1)/pow(temp,2);
            TN /= temp;
            TNN11.col(a-1-risk1_index) = TN;
            TN *= temp;
            risk1_index--;
          }
          else if (survtime(q+1) != survtime(q))
          {
            TN *= H01(risk1_index, 1)/pow(temp,2);
            TNN1.col(a-1-risk1_index) = TN;
            TRNN += TN;
            TN /= H01(risk1_index, 1)/pow(temp,2);
            TN /= temp;
            TNN11.col(a-1-risk1_index) = TN;
            TN *= temp;
            risk1_index--;
          }
          else
          {
            for (q=q+1;q<k;q++)
            {
              temp += exp(MultVV(X2.row(q), gamma1))*FUNEC(0, q);
              for (i=0;i<p1a;i++) N(i) = FUNBEC(i,q);
              TN += exp(MultVV(X2.row(q), gamma1))*N;
              if (q == k-1)
              {
                TN *= H01(risk1_index, 1)/pow(temp,2);
                TNN1.col(a-1-risk1_index) = TN;
                TRNN += TN;
                TN /= H01(risk1_index, 1)/pow(temp,2);
                TN /= temp;
                TNN11.col(a-1-risk1_index) = TN;
                TN *= temp;
                risk1_index--;
                break;
              }
              else if (survtime(q+1) != survtime(q))
              {
                TN *= H01(risk1_index, 1)/pow(temp,2);
                TNN1.col(a-1-risk1_index) = TN;
                TRNN += TN;
                TN /= H01(risk1_index, 1)/pow(temp,2);
                TN /= temp;
                TNN11.col(a-1-risk1_index) = TN;
                TN *= temp;
                risk1_index--;
                break;
              }
              else continue;
            }
          }
    
        }
        else continue;
      }
      TRNN1.col(j) = TRNN;
    }
    else
    {
      if (risk1_index_vtemp>=0)
      {
        if (survtime(j) >= H01(risk1_index_vtemp, 0))
        {
          TRNN1.col(j) = TRNN1.col(j-1);
        }
        else
        {
          risk1_index_vtemp--;
          if (risk1_index_vtemp>=0)
          {
            TRNN = TRNN1.col(j-1);
            TRNN -= TNN1.col(a-1-risk1_index_vtemp-1);
            TRNN1.col(j) = TRNN;
          }
        }
      }
      else
      {
        risk1_index_vtemp=0;
      }
    }
    TRN = TRNN1.col(j);
    
    TRN *= exp(MultVV(X2.row(j), gamma1))*FUNEC(0,j);
    
    if (j==0)
    {
      for (t=0;t<p1a;t++) N(t) = FUNBEC(t,j);
      N *= CumuH01(risk1_index_vttemp)*exp(MultVV(X2.row(j), gamma1));
      TRN -= N;
    }
    else if (survtime(j) >= H01(risk1_index_vttemp,0))
    {
      for (t=0;t<p1a;t++) N(t) = FUNBEC(t,j);
      N *= CumuH01(risk1_index_vttemp)*exp(MultVV(X2.row(j), gamma1));
      TRN -= N;
    }
    else
    {
      risk1_index_vttemp--;
      if (risk1_index_vttemp>=0)
      {
        for (t=0;t<p1a;t++) N(t) = FUNBEC(t,j);
        N *= CumuH01(risk1_index_vttemp)*exp(MultVV(X2.row(j), gamma1));
        TRN -= N;
      }
      else
      {
        risk1_index_vttemp=0;
      }
    }
    
    
    if (survtime(j) >= H01(risk1_index_vtttemp,0))
    {
      if (cmprsk(j) == 1)
      {
        TN = FUNB.col(j) - TNN11.col(a-1-risk1_index_vtttemp);
        TN += TRN;
        for (q=0;q<p1a;q++) S(p1+p1b+p2+q) = TN(q);
      }
      else
      {
        for (q=0;q<p1a;q++) S(p1+p1b+p2+q) = TRN(q);
      }
    }
    else
    {
      risk1_index_vtttemp--;
      if (risk1_index_vtttemp>=0)
      {
        if (cmprsk(j) == 1)
        {
          TN = FUNB.col(j) - TNN11.col(a-1-risk1_index_vtttemp);
          TN += TRN;
          for (q=0;q<p1a;q++) S(p1+p1b+p2+q) = TN(q);
        }
        else
        {
          for (q=0;q<p1a;q++) S(p1+p1b+p2+q) = TRN(q);
        }
      }
      else
      {
        risk1_index_vtttemp=0;
        for (q=0;q<p1a;q++) S(p1+p1b+p2+q) = TRN(q);
      }
    }
    
    
    /* calculate score for vee */
    /*  vee1 */
    if (j == 0)
    {
      temp=0;
    
      TG = 0;
      TRG = 0;
    
      risk1_index=risk1_index_wtemp;
      for (q=j;q<k;q++)
      {
        temp += exp(MultVV(X2.row(q), gamma1))*FUNEC(0, q);
        G = FUNWEC(0,q);
        TG += exp(MultVV(X2.row(q), gamma1))*G;
        if (cmprsk(q) == 1)
        {
          if (q == k-1)
          {
            TG *= H01(risk1_index, 1)/pow(temp,2);
            TGG1(a-1-risk1_index) = TG;
            TRGG += TG;
            TG /= H01(risk1_index, 1)/pow(temp,2);
            TG /= temp;
            TGG11(a-1-risk1_index) = TG;
            TG *= temp;
            risk1_index--;
          }
          else if (survtime(q+1) != survtime(q))
          {
            TG *= H01(risk1_index, 1)/pow(temp,2);
            TGG1(a-1-risk1_index) = TG;
            TRGG += TG;
            TG /= H01(risk1_index, 1)/pow(temp,2);
            TG /= temp;
            TGG11(a-1-risk1_index) = TG;
            TG *= temp;
            risk1_index--;
          }
          else
          {
            for (q=q+1;q<k;q++)
            {
              temp += exp(MultVV(X2.row(q), gamma1))*FUNEC(0, q);
              G = FUNWEC(0,q);
              TG += exp(MultVV(X2.row(q), gamma1))*G;
              if (q == k-1)
              {
                TG *= H01(risk1_index, 1)/pow(temp,2);
                TGG1(a-1-risk1_index) = TG;
                TRGG += TG;
                TG /= H01(risk1_index, 1)/pow(temp,2);
                TG /= temp;
                TGG11(a-1-risk1_index) = TG;
                TG *= temp;
                risk1_index--;
                break;
              }
              else if (survtime(q+1) != survtime(q))
              {
                TG *= H01(risk1_index, 1)/pow(temp,2);
                TGG1(a-1-risk1_index) = TG;
                TRGG += TG;
                TG /= H01(risk1_index, 1)/pow(temp,2);
                TG /= temp;
                TGG11(a-1-risk1_index) = TG;
                TG *= temp;
                risk1_index--;
                break;
              }
              else continue;
            }
          }
    
        }
        else continue;
      }
      TRGG1(j) = TRGG;
    }
    else
    {
      if (risk1_index_wtemp>=0)
      {
        if (survtime(j) >= H01(risk1_index_wtemp, 0))
        {
          TRGG1(j) = TRGG1(j-1);
        }
        else
        {
          risk1_index_wtemp--;
          if (risk1_index_wtemp>=0)
          {
            TRGG = TRGG1(j-1);
            TRGG -= TGG1(a-1-risk1_index_wtemp-1);
            TRGG1(j) = TRGG;
          }
        }
      }
      else
      {
        risk1_index_wtemp=0;
      }
    }
    TRG = TRGG1(j);
    
    TRG *= exp(MultVV(X2.row(j), gamma1))*FUNEC(0,j);
    
    if (j==0)
    {
      G = FUNWEC(0,j);
      G *= CumuH01(risk1_index_wttemp)*exp(MultVV(X2.row(j), gamma1));
      TRG -= G;
    }
    else if (survtime(j) >= H01(risk1_index_wttemp,0))
    {
      G = FUNWEC(0,j);
      G *= CumuH01(risk1_index_wttemp)*exp(MultVV(X2.row(j), gamma1));
      TRG -= G;
    }
    else
    {
      risk1_index_wttemp--;
      if (risk1_index_wttemp>=0)
      {
        G = FUNWEC(0,j);
        G *= CumuH01(risk1_index_wttemp)*exp(MultVV(X2.row(j), gamma1));
        TRG -= G;
      }
      else
      {
        risk1_index_wttemp=0;
      }
    }
    
    
    if (survtime(j) >= H01(risk1_index_wtttemp,0))
    {
      if (cmprsk(j) == 1)
      {
        TG = FUNW(j) - TGG11(a-1-risk1_index_wtttemp);
        TG += TRG;
        S(p1+p1b+p2+p1a) = TG;
      }
      else
      {
        S(p1+p1b+p2+p1a) = TRG;
      }
    }
    else
    {
      risk1_index_wtttemp--;
      if (risk1_index_wtttemp>=0)
      {
        if (cmprsk(j) == 1)
        {
          TG = FUNW(j) - TGG11(a-1-risk1_index_wtttemp);
          TG += TRG;
          S(p1+p1b+p2+p1a) = TG;
        }
        else
        {
          S(p1+p1b+p2+p1a) = TRG;
        }
      }
      else
      {
        risk1_index_wtttemp=0;
        S(p1+p1b+p2+p1a) = TRG;
      }
    }

    /*calculate Sig matrix*/
    
    for(t=0;t<p1a;t++) bsw(t,t) = FUNBS(t,j);
    
    if(p1a>1)
    {
      u=0;
      for(i=1;i<p1a;i++)
      {
        for(t=0;t<p1a-i;t++) {
          bsw(t,i+t) = FUNBS(p1a+u,j);
          bsw(i+t,t) = bsw(t,i+t);
          u++;
        }
      }
    }
    
    for (t=0;t<p1a;t++) {
      bsw(t,p1a) = FUNBW(t,j);
      bsw(p1a,t) = bsw(t,p1a);
    }
    
    bsw(p1a,p1a) = FUNWS(j);
    
    bsw2 = Sig.inverse()*bsw*Sig.inverse() - Sig.inverse();
    
    for (t=0;t<(p1a+1);t++) S(p1+p1b+p2+p1a+1+t) = 0.5*bsw2(t,t);
    u=0;
    for(q=1;q<(p1a+1);q++)
    {
      for(t=0;t<(p1a+1-q);t++) {
        S(p1+p1b+p2+p1a+1+p1a+1+u) = bsw2(t,q+t);
        u++;
      }
    }
  
    SS += MultVVoutprod(S);
  }
  SSinv = SS.inverse();

  Eigen::VectorXd sebeta = Eigen::VectorXd::Zero(p1);
  Eigen::VectorXd setau = Eigen::VectorXd::Zero(p1b);
  Eigen::VectorXd segamma1 = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd sealpha1 = Eigen::VectorXd::Zero(p1a);
  double sevee1=0;
  Eigen::MatrixXd seSig = Eigen::MatrixXd::Zero(p1a+1, p1a+1);


  for (t=0;t<p1;t++) sebeta(t) = sqrt(SSinv(t,t));
  for (t=0;t<p1b;t++) setau(t) = sqrt(SSinv(p1+t,p1+t));
  for (t=0;t<p2;t++) segamma1(t) = sqrt(SSinv(p1+p1b+t,p1+p1b+t));
  for (t=0;t<p1a;t++) sealpha1(t) = sqrt(SSinv(p1+p1b+p2+t,p1+p1b+p2+t));
  sevee1 = sqrt(SSinv(p1+p1b+p2+p1a,p1+p1b+p2+p1a));
  for (t=0;t<(p1a+1);t++) seSig(t,t) = sqrt(SSinv(p1+p1b+p2+p1a+1+t,p1+p1b+p2+p1a+1+t));
  u=0;
  for(q=1;q<(p1a+1);q++)
  {
    for(t=0;t<(p1a+1-q);t++) {
      seSig(t,q+t) = sqrt(SSinv(p1+p1b+p2+p1a+1+p1a+1+u,p1+p1b+p2+p1a+1+p1a+1+u));
      seSig(q+t,t) = seSig(t,q+t);
      u++;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("vcov")=SSinv,
                            Rcpp::Named("seSig")=seSig);
  
  
}
