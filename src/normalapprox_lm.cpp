#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]


Rcpp::List normalApprox_lm(Rcpp::List XList, Rcpp::List YList, Rcpp::List ZList, Eigen::MatrixXd& W,
                           Rcpp::List mdata, Rcpp::List mdataSList,
                           Rcpp::List bList, Eigen::VectorXd sigmaInit, Rcpp::List sigmaiList,
                           Eigen::MatrixXd H01, Eigen::MatrixXd H02, Eigen::VectorXd& survtime, Eigen::VectorXd cmprsk,
                           Eigen::VectorXd& gamma1, Eigen::VectorXd& gamma2, Rcpp::List alphaList,
                           const Eigen::VectorXd& CUH01,
                           const Eigen::VectorXd& CUH02,
                           const Eigen::VectorXd& HAZ01,
                           const Eigen::VectorXd& HAZ02, const Eigen::MatrixXd& Sig,
                           Rcpp::List betaList, Rcpp::List XsList, double s, Eigen::MatrixXd Zs,
                           std::string latAsso){
  
  int numSubj = XList.size();
  int numBio = Rcpp::as<Rcpp::List>(XList[0]).size();
  
  Eigen::VectorXd pVec = Eigen::VectorXd::Zero(numBio);
  int ptotal = 0;
  int p;
  Eigen::VectorXd pREVec = Eigen::VectorXd::Zero(numBio);
  int pREtotal = 0;
  int pRE = 0;
  
  for (int g = 0; g < numBio; g++) {
    p = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(XList[0])[g]).cols();
    pVec(g) = p;
    ptotal += p;
    pRE = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(ZList[0])[g]).cols();
    pREVec(g) = pRE;
    pREtotal += pRE;
  }
  
  // std::cout << "pREVec" << pREVec << std::endl;
  // std::cout << "pREtotal" << pREtotal << std::endl;
  
  
  int index = 0;
  int pREindex = 0;
  Rcpp::List betaNewList;
  Eigen::VectorXd betaFull = Eigen::VectorXd::Zero(ptotal);
  Eigen::VectorXd sigmaVec = Eigen::VectorXd::Zero(numBio);
  
  for(int g = 0; g < numBio; g++){
    
    //~~~~~~~~~~~~
    //
    // BETA
    //
    // ~~~~~~~~~~~~~~~
    
    p = pVec(g);
    pRE= pREVec(g);
    
    Eigen::MatrixXd XVXT = Eigen::MatrixXd::Zero(p, p);
    Eigen::MatrixXd YZBX = Eigen::MatrixXd::Zero(p, 1);
    Eigen::VectorXd betaNew = Eigen::VectorXd::Zero(p);
    
    for (int i = 0; i < numSubj; i++) {
      
      Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
      Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
      Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
      // need to change to matrix for flexible b
      
      Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[g]);
      Eigen::VectorXd Ytemp = Rcpp::as<Eigen::VectorXd>(yListElement[g]);
      Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[g]);
      
      Eigen::VectorXd bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
      Eigen::VectorXd bVecig = bVeci.segment(pREindex, pRE);
      // bVeci contains all random effects for subject i, concatenated across biomarkers.
      // pREindex used to segment it biomarker-by-biomarker.
      
      
      double sigmag = sigmaInit(g);
      
      XVXT = XVXT + Xtemp.transpose() * Xtemp / sigmag; //4x4
      YZBX = YZBX + Xtemp.transpose() * (Ytemp - Ztemp * bVecig) / sigmag; // 4x1
      
    }
    
    betaNew = XVXT.inverse() * YZBX;
    betaFull.segment(index, p) = betaNew;
    // std::cout << "betaFull" <<betaFull << std::endl;
    betaNewList[std::string("beta") + std::to_string(g+1)] = betaNew;
    index += p;
    
    // std::cout << "betaNew " << betaNew << std::endl;
    
    
    //~~~~~~~~~~~~
    //
    // sigma
    //
    // ~~~~~~~~~~~~~~~
    
    double numsig = 0;
    int nijSum = 0;
    
    Eigen::MatrixXd ZZT = Eigen::MatrixXd::Zero(pRE, pRE);
    // pREindex = 0;
    for (int i = 0; i < numSubj; i++) {
      Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
      Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
      Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
      // Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
      
      Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[g]);
      Eigen::VectorXd Ytemp = Rcpp::as<Eigen::VectorXd>(yListElement[g]);
      Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[g]);
      
      Eigen::VectorXd bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
      Eigen::VectorXd bVecig = bVeci.segment(pREindex, pRE);
      
      Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
      Eigen::MatrixXd sigig = sigmai.block(pREindex, pREindex, pRE, pRE);  // sigma i for gth biomarker
      
      Rcpp::List mdataList = Rcpp::as<Rcpp::List>(mdata[g]);
      int numRep = Rcpp::as<int>(mdataList[i]);
      
      for (int nij = 0; nij < numRep; nij++) {
        
        double epsilon = Ytemp(nij) - MultVV(Xtemp.row(nij), betaNew);
        double zb = MultVV(Ztemp.row(nij), bVecig);
        
        ZZT = MultVVoutprod(Ztemp.row(nij));
        Eigen::MatrixXd bbT = MultVVoutprod(bVecig);
        
        // 
        //          std::cout << "zzt" << ZZT << std::endl;
        //         std::cout << "bbt" << bbT << std:: endl;
        
        numsig += pow(epsilon, 2) - 2 * epsilon * zb + (ZZT * (sigig + bbT)).trace();
        
        // if(g == 1){
        // std::cout << "z" << Ztemp.row(nij) << std::endl;
        // std::cout << "bVecig " << bVecig << std::endl;
        // std::cout << "zb " << zb << std::endl;
        // std::cout << "zzt " << ZZT << std::endl;
        // std::cout << "bbt" << bbT << std:: endl;
        //    std::cout << "epsilon" << epsilon << std::endl;
        //    std::cout << "middle part" << 2 * epsilon * zb << std:: endl;
        // std::cout << "trace thing" << (ZZT * (sigig + bbT)).trace() << std::endl;
        // // std::cout << "z" << Ztemp.row(nij) << std::endl;
        // std::cout << "bVecig " << bVecig << std::endl;
        // std::cout << "numsig " <<numsig << std::endl;
        // }
        
      }
      
      nijSum += numRep;
      
    }
    
    sigmaVec(g) = numsig / nijSum;
    
    
    pREindex += pRE;
    // std::cout << "numer" << numsig << std::endl;
    // std::cout << "nij" << nijSum<< std::endl;
    // std::cout << "pREindex" << pREindex << std::endl;
    
  }
  
  
  //~~~~~~~~~~~~
  //
  // SIGMA
  //
  // ~~~~~~~~~~~~~~~
  
  Eigen::MatrixXd SigE = Eigen::MatrixXd::Zero(pREtotal, pREtotal);
  Eigen::MatrixXd numSig = Eigen::MatrixXd::Zero(pREtotal, pREtotal);
  Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(pREtotal);
  
  // int count = 0;
  
  
  for (int i = 0; i < numSubj; i++) {
    
    index = 0;
    
    Eigen::VectorXd bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
    Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    
    numSig += sigmai + MultVVoutprod(bVeci);
    // std::cout << "bVeci" << bVeci << std::endl;
    // std::cout << "sigmai" << sigmai << std::endl;
    // std::cout << "bbT" << MultVVoutprod(bVeci) << std::endl;
    
  }
  
  SigE = numSig / numSubj;
  
  
  
  // ----------------------
  //   Hazard
  // ----------------------
  
  double dem1 = 0;
  double dem2 = 0;
  int a = H01.rows();
  int b = H02.rows();
  int risk1_index = a - 1;
  int risk2_index = b - 1;
  
  
  Eigen::VectorXd alpha1 = Eigen::VectorXd::Zero(numBio);
  Eigen::VectorXd alpha2 = Eigen::VectorXd::Zero(numBio);
  Rcpp::List alphaListElement1 = Rcpp::as<Rcpp::List>(alphaList[0]); // get first  risk
  Rcpp::List alphaListElement2 = Rcpp::as<Rcpp::List>(alphaList[1]); // get second risk
  
  index = 0;
  
  for(int g = 0; g < numBio; g++){
    // pRE = pREVec(g);
    double alpha1g = Rcpp::as<double>(alphaListElement1[g]);
    double alpha2g = Rcpp::as<double>(alphaListElement2[g]);
    
    alpha1[index] = alpha1g;
    alpha2[index] = alpha2g;
    index += 1;
  }
  
  
  index = 0;
  
  // // --- normal ----------------------
  
  Eigen::MatrixXd latent = Eigen::MatrixXd::Zero(numSubj, numBio);
  
  for (int i = 0; i < numSubj; i++) {
    
    Rcpp::List xsListElement = Rcpp::as<Rcpp::List>(XsList[i]);
    Eigen::VectorXd bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
    
    int betaIndex = 0;
    int reIndex = 0;
    
    for (int g = 0; g < numBio; g++) {
      
      int pg = static_cast<int>(pVec(g));
      int qg = static_cast<int>(pREVec(g));
      
// 1 x pg
         // pg x 1
      Eigen::VectorXd b_g = bVeci.segment(reIndex, qg);                   // qg x 1
      
      Eigen::RowVectorXd Zg_s = Zs.block(g, reIndex, 1, qg);              // 1 x qg
      
      
      if (latAsso == "present") {
        Eigen::MatrixXd Xg_s = Rcpp::as<Eigen::MatrixXd>(xsListElement[g]); 
        Eigen::VectorXd beta_g = betaFull.segment(betaIndex, pg);  
        latent(i, g) = (Xg_s * beta_g)(0) + (Zg_s * b_g)(0);
        betaIndex += pg;
      }
      else if (latAsso == "presentlp") {
        latent(i, g) = (Zg_s * b_g)(0);
      }
      
      reIndex += qg;
    }
  }
  
  
  // HAZARD 1
  for (int i = 0; i < numSubj; i++) {
    
    Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    
    // index = 0;
    // int indexX = 0;
    // for(int g = 0; g < numBio; g++){
    //   Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[g]);
    //   Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[g]);
    //   Rcpp::List mdataList = Rcpp::as<Rcpp::List>(mdata[g]);
    //   int numRep = Rcpp::as<int>(mdataList[i]);
    //   
    //   Xs.block(g,indexX, 1, Xtemp.rows()) = Xtemp.row(Xtemp.rows()-1);
    //   Zs.block(g,index, 1,pREVec(g)) = Ztemp.row(Ztemp.rows()-1);
    //   indexX += Xtemp.rows();
    //   index += pREVec(g);
    // }
    
    
    Eigen::VectorXd bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
    
    double muH1, tausq1;
    Eigen::VectorXd latent_i = latent.row(i).transpose();
    muH1 = MultVV(W.row(i), gamma1) + alpha1.dot(latent_i);
    Eigen::MatrixXd BAssociation = Zs * sigmai * Zs.transpose();
    tausq1 = alpha1.transpose() * BAssociation * alpha1;
    dem1 += exp(muH1 + 0.5 * tausq1);
    
    if (cmprsk(i) == 1) {
      //dem += exp(muH1 + 0.5 * tau1);
      
      // check last subject
      if (i == numSubj - 1)
      {
        H01(risk1_index, 2) = H01(risk1_index, 1) / dem1;
        risk1_index--;
      }
      // check if time change
      else if (survtime(i + 1) != survtime(i))
      {
        H01(risk1_index, 2) = H01(risk1_index, 1) / dem1;
        risk1_index--;
      }
      // every other subject
      else
      {
        for (i = i + 1; i < numSubj; i++)
        {
          
          // Xs is wrong
          sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
          
          latent_i = latent.row(i).transpose();
          
          Eigen::VectorXd bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
          
          double muH1, tausq1;
          muH1 = MultVV(W.row(i), gamma1) + alpha1.dot(latent_i);
          Eigen::MatrixXd BAssociation = Zs * sigmai * Zs.transpose();
          
          tausq1 = alpha1.transpose() * BAssociation * alpha1;
          dem1 += exp(muH1 + 0.5 * tausq1);
          
          if (i == numSubj - 1)
          {
            H01(risk1_index, 2) = H01(risk1_index, 1) / dem1;
            risk1_index--;
            break;
          }
          else if (survtime(i + 1) != survtime(i))
          {
            H01(risk1_index, 2) = H01(risk1_index, 1) / dem1;
            risk1_index--;
            break;
          }
          else continue;
        }
      }
      
    }
    else continue;
    
  }
  
  // HAZARD 2
  
  for (int i = 0; i < numSubj; i++) {
    index = 0;
    
    Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
    
    
    Eigen::VectorXd bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
    
    double muH2, tausq2;
    Eigen::VectorXd latent_i = latent.row(i).transpose();
    muH2 = MultVV(W.row(i), gamma2) + alpha2.dot(latent_i);
    Eigen::MatrixXd BAssociation = Zs * sigmai * Zs.transpose();
    
    tausq2 = alpha2.transpose() * BAssociation * alpha2;
    dem2 += exp(muH2 + 0.5 * tausq2);
    
    if (cmprsk(i) == 2) {
      
      // check last subject
      if (i == numSubj - 1)
      {
        H02(risk2_index, 2) = H02(risk2_index, 1) / dem2;
        risk2_index--;
      }
      // check if time change
      else if (survtime(i + 1) != survtime(i))
      {
        H02(risk2_index, 2) = H02(risk2_index, 1) / dem2;
        risk2_index--;
      }
      // every other subject
      else
      {
        for (i = i + 1; i < numSubj; i++)
        {
          
          sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
          
          latent_i = latent.row(i).transpose();;
          
          
          double muH2, tausq2;
          // Eigen::VectorXd latent = Xs * betaFull +Zs * bVeci;
          muH2 = MultVV(W.row(i), gamma2) + alpha2.dot(latent_i);
          Eigen::MatrixXd BAssociation = Zs * sigmai * Zs.transpose();
          
          tausq2 = alpha2.transpose() * BAssociation * alpha2;
          dem2 += exp(muH2 + 0.5 * tausq2);
          
          if (i == numSubj - 1)
          {
            H02(risk2_index, 2) = H02(risk2_index, 1) / dem2;
            risk2_index--;
            break;
          }
          else if (survtime(i + 1) != survtime(i))
          {
            H02(risk2_index, 2) = H02(risk2_index, 1) / dem2;
            risk2_index--;
            break;
          }
          else continue;
        }
      }
    }
    else continue;
  }
  
  /////////////////////////////////////////
  
  /////////////////////////////////////////
  
  // PHI
  
  ////
  
  double scalefH01 = 0;
  double scalefH02 = 0;
  double scalef;
  
  risk1_index = a - 1;
  risk2_index = b - 1;
  
  int dimW = gamma1.size();
  
  Eigen::VectorXd Sw_new = Eigen::VectorXd::Zero(dimW);
  Eigen::VectorXd Sw_inter = Eigen::VectorXd::Zero(dimW);
  Eigen::MatrixXd Sww_new = Eigen::MatrixXd::Zero(dimW, dimW);
  
  Eigen::VectorXd Sl_new = Eigen::VectorXd::Zero(numBio);
  Eigen::VectorXd Sl_inter = Eigen::VectorXd::Zero(numBio);
  Eigen::MatrixXd Sll_new = Eigen::MatrixXd::Zero(numBio, numBio);
  Eigen::MatrixXd Swl_new = Eigen::MatrixXd::Zero(dimW, numBio);
  
  // std::cout<< " dimW " <<  dimW << std::endl;
  
  // 
  Eigen::MatrixXd  wwT = Eigen::MatrixXd::Zero(dimW, dimW);
  Eigen::MatrixXd  llT = Eigen::MatrixXd::Zero(numBio, numBio);
  Eigen::MatrixXd  SwwT = Eigen::MatrixXd::Zero(dimW, dimW);
  Eigen::MatrixXd  SllT = Eigen::MatrixXd::Zero(numBio, numBio);
  Eigen::VectorXd  Sw = Eigen::VectorXd::Zero(dimW);
  Eigen::VectorXd  Sl = Eigen::VectorXd::Zero(numBio);
  Eigen::MatrixXd Swl = Eigen::MatrixXd::Zero(dimW, numBio);
  
  index = 0;
  
  for (int i = 0; i < numSubj; i++) {
    
    Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    
    // index = 0;
    
    Eigen::MatrixXd BAssociation = Zs * sigmai * Zs.transpose(); // numBio x numBio
    Eigen::VectorXd latent_i = latent.row(i).transpose(); //numBio x 1
    double mu1, tausq;
    tausq = alpha1.transpose() * BAssociation * alpha1;
    
    Eigen::VectorXd w = W.row(i);
    Eigen::VectorXd l = BAssociation * alpha1 + latent_i;
    
    Eigen::MatrixXd wl = Eigen::MatrixXd::Zero(dimW, numBio);
    
    wl = w * l.transpose(); //iT
    
    mu1 = MultVV(w, gamma1) + alpha1.dot(latent_i);
    wwT = MultVVoutprod(w);
    llT = MultVVoutprod(BAssociation * alpha1 + latent_i) + BAssociation;
    
    scalef = exp(mu1 + 0.5 * tausq);
    
    // for I
    wwT *= scalef; //exp(mu+tau)wwT
    llT *= scalef; //exp(mu) bbT
    SwwT += wwT; //sum exp(mu)wwT
    SllT += llT; //sum exp(mu)bbT
    wl *= scalef; //exp(mu)bTwT
    Swl += wl; // sum exp(mu)bTwT
    
    // for S
    w *= scalef; //exp(mu)wT
    l *= scalef; //exp(mu)bT
    Sw += w; //sum exp(mu)wT
    Sl += l; //sum exp(mu)
    
    
    if (cmprsk(i) == 1)
    {
      
      if (i == numSubj - 1)
      {
        scalefH01 = H01(risk1_index, 2);
        //SXX *= scalefH01;
        // i
        SwwT *= scalefH01; //haz*exp(mu)wwT
        SllT *= scalefH01; //haz*exp(mu) bbT
        Swl *= scalefH01;
        Sww_new += SwwT;
        Sll_new += SllT;
        Swl_new += Swl;
        
        SwwT /= scalefH01;
        SllT /= scalefH01;
        Swl /= scalefH01;
        
        // s
        Sw *= scalefH01;
        Sl *= scalefH01;
        Sw_new += Sw;
        Sl_new += Sl;
        Sw /= scalefH01;
        Sl /= scalefH01;
        
        
        risk1_index--;
        
      }
      
      else if (survtime(i + 1) != survtime(i))
      {
        scalefH01 = H01(risk1_index, 2);
        
        // i
        SwwT *= scalefH01;
        SllT *= scalefH01;
        Swl *= scalefH01;
        Sww_new += SwwT;
        Sll_new += SllT;
        Swl_new += Swl;
        SwwT /= scalefH01;
        SllT /= scalefH01;
        Swl /= scalefH01;
        
        // s
        Sw *= scalefH01;
        Sl *= scalefH01;
        Sw_new += Sw;
        Sl_new += Sl;
        Sw /= scalefH01;
        Sl /= scalefH01;
        risk1_index--;
      }
      else
      {
        for (i = i + 1; i < numSubj; i++)
        {
          
          Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
          
          
          BAssociation = Zs*sigmai*Zs.transpose();
          latent_i = latent.row(i).transpose();; //numBio x 1
          w = W.row(i);
          l = BAssociation * alpha1 + latent_i;
          
          
          wl = w * l.transpose(); //iT
          mu1 = MultVV(w, gamma1) + alpha1.dot(latent_i);
          wwT =  MultVVoutprod(W.row(i));
          llT =  MultVVoutprod(BAssociation * alpha1 + latent_i) + BAssociation;
          
          tausq = alpha1.transpose() * BAssociation * alpha1;
          scalef = exp(mu1 + 0.5 * tausq);
          
          // for I
          wwT *= scalef; //exp(mu+tau)wwT
          llT *= scalef; //exp(mu) bbT
          SwwT += wwT; //sum exp(mu)wwT
          SllT += llT; //sum exp(mu)bbT
          wl *= scalef; //exp(mu)bTwT
          Swl += wl; // sum exp(mu)bTwT
          
          // for S
          w *= scalef; //exp(mu)wT
          l *= scalef; //exp(mu)bT
          Sw += w; //sum exp(mu)wT
          Sl += l; //sum exp(mu)
          
          
          if (i == numSubj - 1)
          {
            scalefH01 = H01(risk1_index, 2);
            
            // i
            SwwT *= scalefH01;
            SllT *= scalefH01;
            Swl *= scalefH01;
            Sww_new += SwwT;
            Sll_new += SllT;
            Swl_new += Swl;
            SwwT /= scalefH01;
            SllT /= scalefH01;
            Swl /= scalefH01;
            
            // s
            Sw *= scalefH01;
            Sl *= scalefH01;
            Sw_new += Sw;
            Sl_new += Sl;
            Sw /= scalefH01;
            Sl /= scalefH01;
            risk1_index--;
            break;
          }
          else if (survtime(i + 1) != survtime(i))
          {
            scalefH01 = H01(risk1_index, 2);
            
            // i
            SwwT *= scalefH01;
            SllT *= scalefH01;
            Swl *= scalefH01;
            Sww_new += SwwT;
            Sll_new += SllT;
            Swl_new += Swl;
            SwwT /= scalefH01;
            SllT /= scalefH01;
            Swl /= scalefH01;
            
            // s
            Sw *= scalefH01;
            Sl *= scalefH01;
            Sw_new += Sw;
            Sl_new += Sl;
            Sw /= scalefH01;
            Sl /= scalefH01;
            risk1_index--;
            break;
          }
          else continue;
        }
      }
    }
    else continue;
    
    
  }
  
  index = 0;
  
  for (int i = 0; i < numSubj; i++)
  {
    
    
    Eigen::VectorXd latent_i = latent.row(i).transpose();;
    
    if (cmprsk(i) == 1) {
      Sw_inter += W.row(i);
      Sl_inter += latent_i;
    }
  }
  
  
  //NR update
  Eigen::VectorXd Sfull_inter = Eigen::VectorXd::Zero(dimW + numBio);
  Eigen::VectorXd Sfull_new = Eigen::VectorXd::Zero(dimW + numBio);
  Eigen::MatrixXd info = Eigen::MatrixXd::Zero(dimW + numBio, dimW + numBio);
  
  Sfull_inter << Sw_inter, Sl_inter;
  Sfull_new << Sw_new, Sl_new;
  
  // start row, start column, how many rows, how many col
  info.block(0, 0, dimW, dimW) = Sww_new;
  info.block(0, dimW, dimW, numBio) = Swl_new;
  info.block(dimW, 0, numBio, dimW) = Swl_new.transpose();
  info.block(dimW, dimW, numBio, numBio) = Sll_new;
  
  // NR update
  Eigen::VectorXd phi1 = Eigen::VectorXd::Zero(dimW + numBio);
  // std::cout<< "gamma1 " <<  gamma1 << std::endl;
  // std::cout<< "info " << info.inverse() << std::endl;
  phi1 << gamma1, alpha1;
  phi1 += info.inverse() * (Sfull_inter - Sfull_new);
  
  
  // phi 2 ~~
  
  Sw_new = Eigen::VectorXd::Zero(dimW);
  Sw_inter = Eigen::VectorXd::Zero(dimW);
  Sww_new = Eigen::MatrixXd::Zero(dimW, dimW);
  Sl_new = Eigen::VectorXd::Zero(numBio);
  Sl_inter = Eigen::VectorXd::Zero(numBio);
  Sll_new = Eigen::MatrixXd::Zero(numBio, numBio);
  Swl_new = Eigen::MatrixXd::Zero(dimW, numBio);
  
  
  wwT = Eigen::MatrixXd::Zero(dimW, dimW);
  llT = Eigen::MatrixXd::Zero(numBio, numBio);
  SwwT = Eigen::MatrixXd::Zero(dimW, dimW);
  SllT = Eigen::MatrixXd::Zero(numBio, numBio);
  Sw = Eigen::VectorXd::Zero(dimW);
  Sl = Eigen::VectorXd::Zero(numBio);
  Swl = Eigen::MatrixXd::Zero(dimW, numBio);
  
  for (int i = 0; i < numSubj; i++) {
    
    Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    
    Eigen::MatrixXd BAssociation = Zs*sigmai*Zs.transpose(); // numBio x numBio
    Eigen::VectorXd latent_i = latent.row(i).transpose(); //numBio x 1
    double mu2, tausq;
    tausq = alpha2.transpose() * BAssociation * alpha2;
    
    Eigen::VectorXd w = W.row(i);
    Eigen::VectorXd l = BAssociation * alpha2 + latent_i;
    Eigen::MatrixXd wl = Eigen::MatrixXd::Zero(dimW, numBio);
    
    wl = w * l.transpose(); //
    
    mu2 = MultVV(w, gamma2) + alpha2.dot(latent_i);
    wwT = MultVVoutprod(w);
    //llT = MultVVoutprod(l);
    llT = MultVVoutprod(BAssociation * alpha2 + latent_i) + BAssociation;
    
    scalef = exp(mu2 + 0.5 * tausq);
    
    // for I
    wwT *= scalef; //exp(mu+tau)wwT
    llT *= scalef; //exp(mu) bbT
    SwwT += wwT; //sum exp(mu)wwT
    SllT += llT; //sum exp(mu)bbT
    // intersection
    wl *= scalef; //exp(mu)bTwT
    Swl += wl; // sum exp(mu)bTwT
    
    // for S
    w *= scalef; //exp(mu)wT
    l *= scalef; //exp(mu)bT
    Sw += w; //sum exp(mu)wT
    Sl += l; //sum exp(mu)
    
    
    if (cmprsk(i) == 2){
      
      if(i == numSubj - 1){
        
        scalefH02 = H02(risk2_index, 2);
        
        // i
        SwwT *= scalefH02; //haz*exp(mu)wwT
        SllT *= scalefH02; //haz*mm * exp(mu)
        Swl *= scalefH02;
        Sww_new += SwwT;
        Sll_new += SllT;
        Swl_new += Swl;
        
        SwwT /= scalefH02;
        SllT /= scalefH02;
        Swl /= scalefH02;
        
        // s
        Sw *= scalefH02;
        Sl *= scalefH02;
        Sw_new += Sw;
        Sl_new += Sl;
        Sw /= scalefH02;
        Sl /= scalefH02;
        
        risk2_index--;
      }
      
      else if(survtime (i+1) != survtime (i))
      {
        scalefH02 = H02(risk2_index, 2);
        
        // i
        SwwT *= scalefH02; //haz*exp(mu)wwT
        SllT *= scalefH02; //haz*mm * exp(mu)
        Swl *= scalefH02;
        Sww_new += SwwT;
        Sll_new += SllT;
        Swl_new += Swl;
        
        SwwT /= scalefH02;
        SllT /= scalefH02;
        Swl /= scalefH02;
        
        // s
        Sw *= scalefH02;
        Sl *= scalefH02;
        Sw_new += Sw;
        Sl_new += Sl;
        Sw /= scalefH02;
        Sl /= scalefH02;
        
        risk2_index--;
        
      } 
      else{
        for(i = i + 1; i < numSubj; i++){
          Eigen::VectorXd bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
          sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
          
          BAssociation = Zs*sigmai*Zs.transpose();
          latent_i = latent.row(i).transpose();;//numBio x 1
          w = W.row(i);
          l = BAssociation * alpha2 + latent_i;
          
          
          wl = w * l.transpose(); //iT
          mu2 = MultVV(w, gamma2) + alpha2.dot(latent_i);
          wwT =  MultVVoutprod(W.row(i));
          llT = MultVVoutprod(BAssociation * alpha2 + latent_i) + BAssociation;
          
          tausq = alpha2.transpose() * BAssociation * alpha2;
          scalef = exp(mu2 + 0.5 * tausq);
          
          // for I
          wwT *= scalef; //exp(mu+tau)wwT
          llT *= scalef; //exp(mu) bbT
          SwwT += wwT; //sum exp(mu)wwT
          SllT += llT; //sum exp(mu)bbT
          wl *= scalef; //exp(mu)bTwT
          Swl += wl; // sum exp(mu)bTwT
          
          // for S
          w *= scalef; //exp(mu)wT
          l *= scalef; //exp(mu)bT
          Sw += w; //sum exp(mu)wT
          Sl += l; //sum exp(mu)
          
          if(i == numSubj - 1){
            
            scalefH02 = H02(risk2_index, 2);
            
            // i
            SwwT *= scalefH02; //haz*exp(mu)wwT
            SllT *= scalefH02; //haz*mm * exp(mu)
            Swl *= scalefH02;
            Sww_new += SwwT;
            Sll_new += SllT;
            Swl_new += Swl;
            
            SwwT /= scalefH02;
            SllT /= scalefH02;
            Swl /= scalefH02;
            
            // s
            Sw *= scalefH02;
            Sl *= scalefH02;
            Sw_new += Sw;
            Sl_new += Sl;
            Sw /= scalefH02;
            Sl /= scalefH02;
            
            risk2_index--;
            break;
            
          }
          else if (survtime(i + 1) != survtime(i))
          {
            scalefH02 = H02(risk2_index, 2);
            
            // i
            SwwT *= scalefH02; //haz*exp(mu)wwT
            SllT *= scalefH02; //haz*mm * exp(mu)
            Swl *= scalefH02;
            Sww_new += SwwT;
            Sll_new += SllT;
            Swl_new += Swl;
            
            SwwT /= scalefH02;
            SllT /= scalefH02;
            Swl /= scalefH02;
            
            // s
            Sw *= scalefH02;
            Sl *= scalefH02;
            Sw_new += Sw;
            Sl_new += Sl;
            Sw /= scalefH02;
            Sl /= scalefH02;
            
            risk2_index--;
            break;
          }
          else continue;
        }
      }
      
    }
    
    else continue;
  } 
  
  for (int i = 0; i < numSubj; i++){
    Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    
    Eigen::VectorXd latent_i = latent.row(i).transpose();
    if (cmprsk(i) == 2) {
      Sw_inter += W.row(i);
      Sl_inter += latent_i;
      //Sl_inter +=  BAssociation  * alpha1 + bVeci;
    }
  }
  
  
  //NR update
  
  Sfull_inter = Eigen::VectorXd::Zero(dimW + numBio);
  Sfull_new = Eigen::VectorXd::Zero(dimW + numBio);
  info = Eigen::MatrixXd::Zero(dimW + numBio, dimW + numBio);
  
  Sfull_inter << Sw_inter, Sl_inter;
  Sfull_new << Sw_new, Sl_new;
  
  
  // start row, start column, how many rows, how many col
  info.block(0, 0, dimW, dimW) = Sww_new;
  info.block(0, dimW, dimW, numBio) = Swl_new;
  info.block(dimW, 0, numBio, dimW) = Swl_new.transpose();
  info.block(dimW, dimW, numBio, numBio) = Sll_new;
  
  // NR update
  Eigen::VectorXd phi2 = Eigen::VectorXd::Zero(dimW + numBio);
  phi2 << gamma2, alpha2;
  phi2 += info.inverse() * (Sfull_inter - Sfull_new);
  
  // Rcpp::Named("FUNB") = FUNB,
  
  return Rcpp::List::create(Rcpp::Named("beta") = betaFull,
                            Rcpp::Named("betaList") = betaNewList,
                            Rcpp::Named("sigmaVec") = sigmaVec,
                            Rcpp::Named("Sig") = SigE,
                            Rcpp::Named("H01") = H01, Rcpp::Named("H02") = H02,
                            Rcpp::Named("phi1") = phi1, Rcpp::Named("phi2") = phi2);
  
}
