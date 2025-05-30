#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]

Rcpp::List getmvCovSF(const Eigen::VectorXd beta, 
                    const Eigen::VectorXd & gamma1, 
                    const Eigen::VectorXd & alpha1, 
                    const Eigen::MatrixXd & H01, 
                    Rcpp::List sigmaiList,
                    const Eigen::MatrixXd & Sig,
                    const Eigen::VectorXd sigmaVec,
                    Rcpp::List XList, 
                    Rcpp::List YList, 
                    Rcpp::List ZList, 
                    Eigen::MatrixXd& W, 
                    const Eigen::VectorXd & survtime, 
                    const Eigen::VectorXd & cmprsk, 
                    Rcpp::List mdata, Rcpp::List mdataSList, 
                    Rcpp::List bList) {
  
  int numSubj = XList.size();
  int numBio = Rcpp::as<Rcpp::List>(XList[0]).size();

  Eigen::VectorXd pVec = Eigen::VectorXd::Zero(numBio);
  Eigen::VectorXd pREVec = Eigen::VectorXd::Zero(numBio);

  int p, pRE;
  int ptotal = 0;
  int pREtotal = 0;

  for (int g = 0; g < numBio; g++) {
    p = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(XList[0])[g]).cols();
    pVec(g) = p;
    ptotal += p;
    pRE = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(ZList[0])[g]).cols();
    pREVec(g) = pRE;
    pREtotal += pRE;
  }

  // 
  int d = beta.size() + gamma1.size()  + alpha1.size() + Sig.cols()*(Sig.cols() + 1)/2 + numBio;
  int a = H01.rows();

  int i,q,j,t,u;

  double temp,temp1;

  double tausq;//, FUNBEC;

  int i2;


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

  Eigen::VectorXd CumuH01 = Eigen::VectorXd::Zero(a);

  temp1=0;
  for (i=0;i<a;i++) {
    temp1 += H01(i, 2);
    CumuH01(i) = temp1;
  }
  temp1=0;

  double epsilon=0;



  double qqsigma=0;



  // int p1 = beta.size();
  // int p2 = gamma1.size();
  //
  // biomarker loop here
  // 
  // 
  Eigen::VectorXd SZ = Eigen::VectorXd::Zero(p);
  // Eigen::VectorXd SZ1 = Eigen::VectorXd::Zero(p1);
  // Eigen::MatrixXd SZZ = Eigen::MatrixXd::Zero(p1, pREtotal);

  Eigen::VectorXd X = Eigen::VectorXd::Zero(gamma1.size());
  Eigen::VectorXd SX = Eigen::VectorXd::Zero(gamma1.size());
  Eigen::VectorXd SX1 = Eigen::VectorXd::Zero(gamma1.size());
  Eigen::VectorXd SRX = Eigen::VectorXd::Zero(gamma1.size());
  Eigen::VectorXd SRXX = Eigen::VectorXd::Zero(gamma1.size());
  Eigen::MatrixXd SXX1 = Eigen::MatrixXd::Zero(gamma1.size(), a);
  Eigen::MatrixXd SXX11 = Eigen::MatrixXd::Zero(gamma1.size(), a);
  Eigen::MatrixXd SRXX1 = Eigen::MatrixXd::Zero(gamma1.size(), numSubj);
  // 
  // 
  // // NEED TO ADAPT FOR DIFFERENT NUMBER OF RANDOM EFFECTS
  Eigen::VectorXd N = Eigen::VectorXd::Zero(pREtotal);
  Eigen::VectorXd TN = Eigen::VectorXd::Zero(pREtotal);
  Eigen::VectorXd TN1 = Eigen::VectorXd::Zero(pREtotal);
  Eigen::VectorXd TRN = Eigen::VectorXd::Zero(pREtotal);
  Eigen::VectorXd TRNN = Eigen::VectorXd::Zero(pREtotal);
  Eigen::MatrixXd TRNN1 = Eigen::MatrixXd::Zero(pREtotal,numSubj);
  Eigen::MatrixXd TRNN2 = Eigen::MatrixXd::Zero(pREtotal,numSubj);
  Eigen::MatrixXd TNN1 = Eigen::MatrixXd::Zero(pREtotal,a);
  Eigen::MatrixXd TNN11 = Eigen::MatrixXd::Zero(pREtotal,a);

  Eigen::MatrixXd bs = Eigen::MatrixXd::Zero(pREtotal,pREtotal);
  Eigen::MatrixXd ZZT = Eigen::MatrixXd::Zero(pRE,pRE);
  Eigen::MatrixXd bs2 = Eigen::MatrixXd::Zero(pREtotal,pREtotal);


  double r;//, epsilon;
  //int numSubj = XList.size();
  //int numBio = Rcpp::as<Rcpp::List>(XList[0]).size();

  // TEMPORARY MAKE SURE TO COMMENT OUT LATER!!!
  // 
  // 
  // 
  // // change to multiple biomarker
  // 
  int index = 0;
  int pREindex = 0;
  int agIndex = 0;
  double num;
  int SEindex = 0;

  //    Eigen::MatrixXd FUNEC = Eigen::MatrixXd::Zero(2,numSubj);
  //
  // //    //b
  Eigen::MatrixXd FUNB = Eigen::MatrixXd::Zero(pREtotal,numSubj);
  Eigen::MatrixXd FUNBS = Eigen::MatrixXd::Zero(pREtotal*(pREtotal+1)/2,numSubj); //bbt
  // Eigen::VectorXd bifull = Eigen::VectorXd::Zero(pREtotal);
  Eigen::MatrixXd bVeci  = Eigen::VectorXd::Zero(pREtotal);
  Eigen::MatrixXd BAssociation = Eigen::MatrixXd::Zero(pREtotal,pREtotal);
  Eigen::MatrixXd FUNEC = Eigen::MatrixXd::Zero(2, numSubj); //exp(alpha*b)
  //     //bexp(alpha*b)
  // mult pRetotal by num CR
  Eigen::MatrixXd FUNBEC = Eigen::MatrixXd::Zero(pREtotal,numSubj); // 2 for number of competing risks
  // 
  // 
  // //     //bbTexp(alpha*b)
      // Eigen::MatrixXd FUNBSEC = Eigen::MatrixXd::Zero(2*pRE*(pRE+1)/2,k);

  Eigen::VectorXd FUNBECvec = Eigen::VectorXd::Zero(pREtotal);
  // 
  // 
  for(int i = 0; i < numSubj; i++){

    bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
    pREindex = 0;

    // FUNB.col(i) = bifull;
    FUNB.col(i) = bVeci;

    // bVeci = bifull;
    // latent
    Eigen::VectorXd latent = bVeci;
    BAssociation = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]); // this part need ot make slexible
    Eigen::MatrixXd sigi = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);

    Eigen::MatrixXd FUNBSmat =  sigi + MultVVoutprod(bVeci); // pREtotalxpRetotal
    Eigen::VectorXd FUNBSvec =  Eigen::VectorXd::Zero(pREtotal*(pREtotal+1)/2-FUNBSmat.diagonal().size()); // 6

    int FUNBSind= 0;

    for (int row = 1; row < pREtotal; row++) {
      for (int col = 0; col < pREtotal-row; col++) {
        FUNBSvec(FUNBSind) = FUNBSmat(col, col+row);
        FUNBSind++;
      }
    }

    FUNBS.col(i) << FUNBSmat.diagonal(), FUNBSvec;

    tausq = alpha1.transpose() * BAssociation * alpha1;
    FUNEC(0,i) = exp(MultVV(alpha1, latent) + 0.5 * tausq);
    FUNBEC.col(i).segment(index, pREtotal) = exp(MultVV(alpha1, latent) + 0.5 * tausq) * (BAssociation * alpha1 + latent);
    // 
  }

  SZ = Eigen::VectorXd::Zero(p);

  Eigen::VectorXd betaVec = Eigen::VectorXd::Zero(ptotal);
  Eigen::VectorXd qqsigmaVec = Eigen::VectorXd::Zero(numBio);


  for (int i = 0; i < numSubj; i++) {

    index = 0;
    S = Eigen::VectorXd::Zero(d); // vector that stores observed score vector of each subject; when looping, will clear up for each subject!!


    // ~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // BETA
    //
    // ~~~~~~~~~~~~~~~~~~~~~~~~


    int betaIndex = 0;
    //int sigIndex = 0;
    int bIndex = 0;
    for(int g = 0; g < numBio; g++){
      double sigma = sigmaVec(g);
      pRE = pREVec(g);
      p = pVec(g);
      num = 0;
      Rcpp::List mdataList = Rcpp::as<Rcpp::List>(mdata[g]);
      
      
      int numRep = Rcpp::as<int>(mdataList[i]);
      SZ = Eigen::VectorXd::Zero(p);

      Eigen::VectorXd betag = beta.segment(betaIndex,p);
      Eigen::VectorXd bVecig = FUNB.col(i).segment(bIndex,pRE);

      Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
      Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
      Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);


      Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[g]);
      Eigen::VectorXd Ytemp = Rcpp::as<Eigen::VectorXd>(yListElement[g]); // need to change to vector later
      Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[g]);


      for (int nij = 0; nij < numRep; nij++) {


        // need to generalize this
        double r = Ytemp(nij) - MultVV(Xtemp.row(nij), betag); //here
        double zb = MultVV(Ztemp.row(nij), bVecig); // here
    
        epsilon = r - zb;
        // beta calc
        SZ += epsilon*Xtemp.row(nij);

        // sigma calc
        Eigen::MatrixXd ZZT = MultVVoutprod(Ztemp.row(nij));
        Eigen::MatrixXd bbT = MultVVoutprod(bVecig);

        Eigen::MatrixXd sigig = BAssociation.block(bIndex, bIndex, pRE, pRE); // THIS BLOCK NEEDS TO make sure it's subject specific
        num += pow(r, 2) - 2 * r * zb + (ZZT * (sigig + bbT)).trace();

      }

      SZ /= sigma;
      qqsigma = num/(2*pow(sigma,2)) - numRep/(2*sigma);

      betaVec.segment(betaIndex,p) = SZ;
      qqsigmaVec(g) = qqsigma;
      betaIndex += p;
      bIndex += pRE;
    }


    S.segment(index,ptotal + numBio) << betaVec, qqsigmaVec;
  
    index += ptotal + numBio;

    ///////////////////////////
    
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // GAMMA
    //
    // ~~~~~~~~~~~~~~~~~~~~~~~~
    
    /* calculate score for gamma */
    if (i == 0){
      temp=0;
      risk1_index=risk1_index_temp;
      for (i2=i;i2<numSubj;i2++) // double check this index
      {

        Eigen::VectorXd w = W.row(i2);
        temp += exp(MultVV(w, gamma1))*FUNEC(0,i2);
        SX += exp(MultVV(w,gamma1))*FUNEC(0,i2)*w;

        if (cmprsk(i2) == 1)
        {
          if (i2 == numSubj-1)
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
          else if (survtime(i2+1) != survtime(i2))
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
            for (i2=i2+1;i2<numSubj;i2++) // subject number
            {

              Eigen::VectorXd w = W.row(i2);
              temp += exp(MultVV(w, gamma1))*FUNEC(0,i2);
              SX += exp(MultVV(w,gamma1))*FUNEC*w;


              if (i2 == numSubj-1)
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
              else if (survtime(i2+1) != survtime(i2))
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
      SRXX1.col(i) = SRXX;
    }
    else
    {
      if (risk1_index_temp>=0)
      {
        if (survtime(i) >= H01(risk1_index_temp, 0))
        {
          SRXX1.col(i) = SRXX1.col(i-1);
        }
        else
        {
          risk1_index_temp--;
          if (risk1_index_temp>=0)
          {
            SRXX = SRXX1.col(i-1);
            SRXX -= SXX1.col(a-1-risk1_index_temp-1);
            SRXX1.col(i) = SRXX;
          }
        }
      }
      else
      {
        risk1_index_temp=0;
      }
    }
    SRX = SRXX1.col(i);

    if (i==0)
    {
      SRX -= CumuH01(risk1_index_ttemp)*W.row(i);
    }
    else if (survtime(i) >= H01(risk1_index_ttemp, 0))
    {
      SRX -= CumuH01(risk1_index_ttemp)*W.row(i);
    }
    else
    {
      risk1_index_ttemp--;
      if (risk1_index_ttemp>=0)
      {
        SRX -= CumuH01(risk1_index_ttemp)*W.row(i);
      }
      else
      {
        SRX = SRXX1.col(i);
        risk1_index_ttemp=0;
      }
    }


    SRX*= exp(MultVV(W.row(i), gamma1))*FUNEC(0,i);

    if (survtime(i) >= H01(risk1_index_tttemp, 0))
    {
      if (cmprsk(i) == 1)
      {
        X = W.row(i);
        X -= SXX11.col(a-1-risk1_index_tttemp);
        X += SRX;
        S.segment(index, gamma1.size()) = X;
      }
      else
      {
        S.segment(index, gamma1.size()) = SRX;
      }
    }
    else
    {
      risk1_index_tttemp--;
      if (risk1_index_tttemp>=0)
      {
        if (cmprsk(i) == 1)
        {
          X = W.row(i);
          X -= SXX11.col(a-1-risk1_index_tttemp);
          X += SRX;
          S.segment(index, gamma1.size()) = X;
        }
        else
        {
          S.segment(index, gamma1.size()) = SRX;
        }
      }
      else
      {
        risk1_index_tttemp=0;
        S.segment(index, gamma1.size()) = SRX;
      }
    }

    // std::cout << "gamma1 " << S.segment(index, gamma1.size())<< std::endl;
    // std::cout << "alpha1full " << S << std::endl;
    index += gamma1.size();



    // ~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // ALPHA
    //
    // ~~~~~~~~~~~~~~~~~~~~~~~~


    /* calculate score for alpha */
    /*  alpha1 */
    if (i == 0){
      temp=0;

      TN = Eigen::VectorXd::Zero(pREtotal);
      TRN = Eigen::VectorXd::Zero(pREtotal);

      risk1_index=risk1_index_vtemp;
      for (i2=i;i2<numSubj;i2++){
        Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(pREtotal);
        Eigen::MatrixXd BAssociation = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i2]);
        Eigen::VectorXd latent = bVeci;

        tausq = alpha1.transpose() * BAssociation * alpha1;

        temp += exp(MultVV(W.row(i2), gamma1))*FUNEC(0, i2);

        N = FUNBEC.col(i2).segment(agIndex, pREtotal);


        TN += exp(MultVV(W.row(i2), gamma1))*N; // E(bTexp(aTb))xexp(WTgamma)


        if (cmprsk(i2) == 1)
        {
          if (i2 == numSubj-1)
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
          else if (survtime(i2+1) != survtime(i2))
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
            for (i2=i2+1;i2<numSubj;i2++)
            {
              temp += exp(MultVV(W.row(i2), gamma1))*FUNEC(0, i2);
              for (int ind=0;ind<pREtotal;ind++) N(ind) = FUNBEC(ind,i2);
              TN += exp(MultVV(W.row(i2), gamma1))*N;
              if (i2 == numSubj-1)
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
              else if (survtime(i2+1) != survtime(i2))
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
      TRNN1.col(i) = TRNN;
    }
    else{
      if (risk1_index_vtemp>=0)
      {
        if (survtime(i) >= H01(risk1_index_vtemp, 0))
        {
          TRNN1.col(i) = TRNN1.col(i-1);
        }
        else
        {
          risk1_index_vtemp--;
          if (risk1_index_vtemp>=0)
          {
            TRNN = TRNN1.col(i-1);
            TRNN -= TNN1.col(a-1-risk1_index_vtemp-1);
            TRNN1.col(i) = TRNN;
          }
        }
      }
      else
      {
        risk1_index_vtemp=0;
      }
    }
    TRN = TRNN1.col(i);

    TRN *= exp(MultVV(W.row(i), gamma1))*FUNEC(0,i);

    if (i==0)
    {
      for (t=0;t<pREtotal;t++) N(t) = FUNBEC(t,i);
      N *= CumuH01(risk1_index_vttemp)*exp(MultVV(W.row(i), gamma1));
      TRN -= N;
    }
    else if (survtime(i) >= H01(risk1_index_vttemp,0))
    {
      for (t=0;t<pREtotal;t++) N(t) = FUNBEC(t,i);
      N *= CumuH01(risk1_index_vttemp)*exp(MultVV(W.row(i), gamma1));
      TRN -= N;
    }
    else
    {
      risk1_index_vttemp--;
      if (risk1_index_vttemp>=0)
      {
        for (t=0;t<pREtotal;t++) N(t) = FUNBEC(t,i);
        N *= CumuH01(risk1_index_vttemp)*exp(MultVV(W.row(i), gamma1));
        TRN -= N;
      }
      else
      {
        risk1_index_vttemp=0;
      }
    }


    if (survtime(i) >= H01(risk1_index_vtttemp,0))
    {
      if (cmprsk(i) == 1)
      {
        TN = FUNB.col(i) - TNN11.col(a-1-risk1_index_vtttemp);
        TN += TRN;

        for (q=0;q<pREtotal;q++) S(index + q) = TN(q);

      }
      else
      {
        for (q=0;q<pREtotal;q++) S(index + q) = TRN(q);
      }
    }
    else
    {

      risk1_index_vtttemp--;
      if (risk1_index_vtttemp>=0)
      {
        if (cmprsk(i) == 1)
        {
          TN = FUNB.col(i) - TNN11.col(a-1-risk1_index_vtttemp);
          TN += TRN;
          for (q=0;q<pREtotal;q++) S(index + q) = TN(q);
        }
        else
        {
          for (q=0;q<pREtotal;q++) S(index + q) = TRN(q);
        }
      }
      else
      {
        risk1_index_vtttemp=0;
        for (q=0;q<pREtotal;q++) S(index + q) = TRN(q);
      }
    }

    // std::cout << "alpha1 " << S(index + q) << std::endl;
    // std::cout << "alpha1full " << S << std::endl;
    index += pREtotal; // double check this

    // std::cout << "alpha2 " << S(index + q) << std::endl;
    // std::cout << "alpha2full " << S << std::endl;


    // agIndex+= numBio;
    // pREindex +=pRE;
 



    // ~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // SIGMA
    //
    // ~~~~~~~~~~~~~~~~~~~~~~~~



    for(t=0;t<pREtotal;t++) bs(t,t) = FUNBS(t,i);// needs to therefore be in loop

    if(pREtotal>1)
    {
      u=0;
      for(int ind=1;ind<pREtotal;ind++)
      {
        for(t=0;t<pREtotal-ind;t++) {
          bs(t,ind+t) = FUNBS(pREtotal+u,i);
          bs(ind + t,t) = bs(t,ind+t);
          u++;
        }
      }
    }


    bs = Sig.inverse()*bs*Sig.inverse() - Sig.inverse();

    for (t=0;t<pREtotal;t++) S(index+t) = 0.5*bs(t,t);
    index += pREtotal;

    u=0;
    for(q=1;q<pREtotal;q++)
    {
      for(t=0;t<(pREtotal-q);t++) {
        S(index+u) = bs(t, q+t);
        u++;
      }
    }


    SS += MultVVoutprod(S);



  }

  SSinv = SS.inverse();


  Eigen::VectorXd sebeta = Eigen::VectorXd::Zero(ptotal);
  Eigen::VectorXd segamma1 = Eigen::VectorXd::Zero(gamma1.size());
  Eigen::VectorXd sealpha1 = Eigen::VectorXd::Zero(pREtotal);
  Eigen::VectorXd sesigma = Eigen::VectorXd::Zero(numBio);
  Eigen::MatrixXd seSig = Eigen::MatrixXd::Zero(pREtotal, pREtotal);



  index = 0;

  // beta
  for (t=0;t<ptotal;t++) sebeta(t) = sqrt(SSinv(t,t));
  index += ptotal;

  // sigma
  for (t=0;t< numBio;t++)  sesigma(t) = sqrt(SSinv(index + t,index + t));
  index += numBio;

  // gamma
  for (t=0;t<gamma1.size();t++) segamma1(t) = sqrt(SSinv(index+t,index+t));
  index += gamma1.size();

  // alpha
  for (t=0;t<pREtotal;t++) sealpha1(t) = sqrt(SSinv(index+t,index+t));
  index += pREtotal;

  // sigma

  for (t=0;t<pREtotal;t++){ seSig(t,t) = sqrt(SSinv(index+t,index+t));
  }
  index += pREtotal;
  u=0;
  for(q=1;q<pREtotal;q++)
  {
    for(t=0;t<(pREtotal-q);t++) {
      seSig(t,q+t) = sqrt(SSinv(index+u,index+u));
      seSig(q+t,t) = seSig(t,q+t);
      u++;
    }
  }


  return Rcpp::List::create(Rcpp::Named("vcov")=SSinv,
                            Rcpp::Named("sebeta")=sebeta,
                            Rcpp::Named("sesigma")=sesigma,
                            Rcpp::Named("segamma1")=segamma1,
                            Rcpp::Named("sealpha1")=sealpha1,
                            Rcpp::Named("seSig")=seSig,
                            Rcpp::Named("FUNBS")= FUNBS,
                            Rcpp::Named("FUNBEC")= FUNBEC,
                            Rcpp::Named("S")=SS);
  
}
