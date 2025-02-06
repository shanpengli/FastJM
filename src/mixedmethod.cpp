#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]


Rcpp::List getQuadMix(Rcpp::List XList, Rcpp::List YList, Rcpp::List ZList, Eigen::MatrixXd& W,
                      Rcpp::List mdata, Rcpp::List mdataSList, 
                      Rcpp::List bList, Rcpp::List sigmaList, Rcpp::List sigmaiList,
                      Eigen::VectorXd weight, Eigen::VectorXd absc,
                      Eigen::MatrixXd H01, Eigen::MatrixXd H02, Eigen::VectorXd& survtime, Eigen::VectorXd cmprsk,
                      Eigen::VectorXd& gamma1, Eigen::VectorXd& gamma2, Rcpp::List alphaList,
                      const Eigen::VectorXd& CUH01,
                      const Eigen::VectorXd& CUH02,
                      const Eigen::VectorXd& HAZ01,
                      const Eigen::VectorXd& HAZ02, const Eigen::MatrixXd& Sig,
                      Rcpp::List betaList){
  
  
  // Eigen::MatrixXd H01q = H01;
  // Eigen::MatrixXd H02q = H02;
  int numWeight = weight.size();
  
  
  
  int numSubj = XList.size();
  int numBio = Rcpp::as<Rcpp::List>(XList[0]).size();
  int p1 = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(XList[0])[0]).cols();
  int p2 = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(XList[0])[1]).cols();
  Eigen::MatrixXd XVXT = Eigen::MatrixXd::Zero(p1, p1);
  Eigen::MatrixXd YZBX = Eigen::MatrixXd::Zero(p1, 1);
  Eigen::VectorXd beta1New = Eigen::VectorXd::Zero(p1);
  Eigen::VectorXd beta2New = Eigen::VectorXd::Zero(p2);
  Eigen::VectorXd betaNew = Eigen::VectorXd::Zero(p1 + p2);
  
  /*******************
   *  ADAPT-GAUSS BETA CALCULATION
   ********************/
  //Rcpp::List mdataList = Rcpp::as<Rcpp::List>(mdata[0]);
  ////int numSubj = mdataList.size();
  //int p1a = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(ZList[0])[0]).cols();
  //int p2a = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(ZList[0])[1]).cols();
  
  //double dem, cuh01, cuh02, haz01, haz02, xgamma1, xgamma2, temp, mu, zb;
  //Eigen::VectorXd c(p1a + p2a);
  //Eigen::VectorXd bi(p1a + p2a);
  //Eigen::VectorXd bi2(p1a + p2a);
  //Eigen::VectorXd weightbi(p1a + p2a);
  //Eigen::VectorXd ri(p1a + p2a);
  //Eigen::VectorXd rii(p1a + p2a);
  //Eigen::MatrixXd Hi(p1a + p2a, p1a + p2a);
  //Eigen::MatrixXd Hi2(p1a + p2a, p1a + p2a);
  
  //Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig.inverse(), Eigen::ComputeThinU | Eigen::ComputeThinV);
  //Eigen::VectorXd eigenSQ = svd.singularValues();
  ////int i, j, q, t, db;
  //for (int i = 0; i < eigenSQ.size(); i++) {
  //	eigenSQ(i) = sqrt(eigenSQ(i));
  //}
  //Eigen::MatrixXd SigSQRT = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  ////std::cout << "Sigsqrtr\n " << SigSQRT << std::endl;
  //Eigen::MatrixXd FUNB = Eigen::MatrixXd::Zero(p1a+p2a, numSubj);
  
  //int point = wsmatrix.rows();
  //int numRep;
  //for (int i = 0; i < numSubj; i++)
  //{
  //	dem = 0;
  //	numRep = mdataList(i);
  //	cuh01 = CUH01(i);
  //	cuh02 = CUH02(i);
  //	haz01 = HAZ01(i);
  //	haz02 = HAZ02(i);
  
  //	xgamma1 = MultVV(W.row(i), gamma1);
  //	xgamma2 = MultVV(W.row(i), gamma2);
  
  //	Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
  
  //	Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(p1a + p2a);
  
  //	//Eigen::MatrixXd bVeci = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
  //	Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
  //	//Eigen::MatrixXd sigig = sigmai.block(0, 0, p1a, p1a);
  
  //	Eigen::VectorXd bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
  //	Eigen::VectorXd bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
  
  //
  //	bVeci << bVec1, bVec2;
  
  //	for (int ct = 0; ct < p1a+p2a; ct++) Hi.row(ct) = sigmai.row(ct);
  //       //Hi.row(ct) = Poscov.row(j * p1a + ct);
  //		
  //	//std::cout << "Hi\n " << Hi << std::endl;
  
  
  //	Eigen::JacobiSVD<Eigen::MatrixXd> svd(Hi, Eigen::ComputeThinU | Eigen::ComputeThinV);
  //	Eigen::VectorXd eigenSQ = svd.singularValues();
  //	for (int ct = 0; ct < eigenSQ.size(); ct++) {
  //		eigenSQ(ct) = sqrt(eigenSQ(ct));
  //	}
  
  //	Hi2 = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  //	//std::cout << "Hi2\n " << Hi2 << std::endl;
  
  //	for (int db = 0; db < point; db++) {
  //		c = xsmatrix.row(db);
  //		weightbi = wsmatrix.row(db);
  //		////bii = Posbi.row(j);
  //		 //4x1 4x4
  //		ri = bVeci + sqrt(2) * Hi2 * c;
  //		rii = SigSQRT * ri;
  //		//temp = exp(705);
  //           temp = 0;
  //           
  //		//biomarker; take first 2 values of ri. line 140-end of loop;; take first 2 ri, bottom 2 ri for each biomarker
  //		// fopr survval/prior dist of re do whole thing
  
  //		for (int ct = 0; ct < p1a+p2a; ct++) c(ct) = ri(ct);
  
  //		//std::cout << "c\n" << c << std::endl;
  //		//std::cout << "c\n" << ri << std::endl;
  
  //		Rcpp::List mdataList = Rcpp::as<Rcpp::List>(mdata[0]);
  //		int numRep = Rcpp::as<int>(mdataList[i]);
  //		Rcpp::List mdataS = Rcpp::as<Rcpp::List>(mdataSList[0]);
  
  //		
  
  //		for (int g = 0; g < numBio; g++) {
  
  //			Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
  //			Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
  //			Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
  //			//Rcpp::List betag = Rcpp::as<Rcpp::List>(BetaList[g]);
  //			Eigen::VectorXd betag = Rcpp::as<Eigen::VectorXd>(betaList[g]);
  //			Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[g]);
  //			Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[g]);
  //			Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[g]);
  
  //			double sigma = Rcpp::as<double>(sigmaList[g]);
  
  //			int index = 0;
  //			Eigen::VectorXd cbio = c.segment(index, Ztemp.cols());
  //			index += Ztemp.cols();
  //			
  //			for (int nij = 0; nij < numRep; nij++) {
  //				mu = MultVV(Xtemp.row(nij), betag);
  //				//std::cout << "xtemp\n" << Xtemp.row(nij) << std::endl;
  //				//std::cout << "beta" << beta1 << std::endl;
  
  //				//std::cout << "ind" << index << std::endl;
  //				//std::cout << "b" << bi << std::endl;
  
  //				zb = MultVV(Ztemp.row(nij), cbio);
  //				//double tempjr = temp;
  //				//temp *= exp(-1 / (2 * sigma) * pow((Ytemp(nij) - mu - zb), 2));
  //                   temp += -1 / (2 * sigma) * pow((Ytemp(nij) - mu - zb), 2);
  //           
  
  //			}
  
  //		}
  
  //		Rcpp::List alphaListElement1 = Rcpp::as<Rcpp::List>(alphaList[0]); // get first  risk
  //		Rcpp::List alphaListElement2 = Rcpp::as<Rcpp::List>(alphaList[1]); // get second risk
  //		Eigen::VectorXd alpha11 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[0]); // get alpha1// risk 1 b1
  //		Eigen::VectorXd alpha12 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[1]); // get alpha2// risk 1 b2
  //		Eigen::VectorXd alpha21 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[0]); // get alpha1// risk 2 b1
  //		Eigen::VectorXd alpha22 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[1]); // get alpha2// risk 2 b2
  //		//Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
  
  //		Eigen::VectorXd alpha1 = Eigen::VectorXd::Zero(p1a + p2a);
  //		Eigen::VectorXd alpha2 = Eigen::VectorXd::Zero(p1a + p2a);
  
  //		alpha1 << alpha11, alpha12; // diff bio, same alpha
  //		alpha2 << alpha21, alpha22; // diff bio, same alpha
  //		
  
  //		//if (cmprsk(i) == 1)  temp *= haz01 * exp(xgamma1 + MultVV(alpha1, bi));
  //		//if (cmprsk(i) == 2)  temp *= haz02 * exp(xgamma2 + MultVV(alpha2, bi));
  //		if (cmprsk(i) == 1) { 
  //			temp += log(haz01) + (xgamma1 + MultVV(alpha1, c)); 
  //			//std::cout << "log" << log(haz01) + (xgamma1 + MultVV(alpha1, bi)) << std::endl;
  //		}
  //		if (cmprsk(i) == 2) {
  //			temp += log(haz02) + (xgamma2 + MultVV(alpha2, c));
  //			//std::cout << "log" << log(haz02) + (xgamma2 + MultVV(alpha2, bi)) << std::endl;
  //		}
  //		
  //		//std::cout << "temp2\n" << temp << std::endl;
  
  //		//temp *= exp(0 - cuh01 * exp(xgamma1 + MultVV(alpha1, bi)) - cuh02 * exp(xgamma2 + MultVV(alpha2, bi)));
  //           temp += - cuh01 * exp(xgamma1 + MultVV(alpha1, c)) - cuh02 * exp(xgamma2 + MultVV(alpha2, c));
  
  //		//std::cout << "temp3\n" << temp << std::endl;
  
  //		for (int ct = 0; ct < p1a+p2a; ct++) temp *= weightbi(ct);
  //		// prior distribution of RE; evaluvate whole vector
  //		bi2 = xsmatrix.row(db);
  //		//temp *= exp(-pow(rii.norm(), 2) / 2) * exp(pow(bi2.norm(), 2));
  //           temp += (-pow(rii.norm(), 2) / 2) + pow(bi2.norm(), 2);
  //                                                   
  //           temp = exp(temp);
  //           
  //		dem += temp;
  //		//calculate h(bi)
  //		FUNB.col(i) += temp * c;
  
  //		
  
  //	}                                                            
  //                                 
  //       FUNB.col(i)/=dem;
  //}
  //
  //for (int i = 0; i < numSubj; i++) {
  //	Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
  //	Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
  //	Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
  //	Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
  
  //	Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[0]);
  //	Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[0]);
  //	Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[0]);
  
  //	double sigma = Rcpp::as<double>(sigmaList[0]);
  
  //	Eigen::MatrixXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
  
  //	XVXT = XVXT + Xtemp.transpose() * Xtemp / sigma; //4x4
  
  //	YZBX = YZBX + Xtemp.transpose() * (Ytemp - Ztemp * FUNB.col(i)) / sigma; // 4x1
  
  //}
  
  
  
  //beta1New = XVXT.inverse() * YZBX;
  
  //XVXT = Eigen::MatrixXd::Zero(p1, p1);
  //YZBX = Eigen::MatrixXd::Zero(p1, 1);
  
  //for (int i = 0; i < numSubj; i++) {
  
  //	Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
  //	Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
  //	Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
  //	Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
  
  
  //	Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[1]);
  //	Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[1]);
  //	Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[1]);
  
  //	double sigma = Rcpp::as<double>(sigmaList[1]);
  //	Eigen::MatrixXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
  
  //	XVXT = XVXT + Xtemp.transpose() * Xtemp / sigma; //4x4
  //	YZBX = YZBX + Xtemp.transpose() * (Ytemp - Ztemp * FUNB.col(i)) / sigma; // 4x1
  
  
  
  //}
  
  //beta2New = XVXT.inverse() * YZBX;
  
  //betaNew << beta1New, beta2New;
  
  /*******************
   *  BETA(NA) CALCULATION
   ********************/
  
  
  
  for (int i = 0; i < numSubj; i++) {
    // Ensure that the current subject exists in XList, YList, ZList
    //if (i < XList.size() && i < YList.size() && i < ZList.size()) {
    Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
    Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
    Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
    Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
    
    Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[0]);
    Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[0]);
    Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[0]);
    
    double sigma = Rcpp::as<double>(sigmaList[0]);
    
    Eigen::MatrixXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
    
    XVXT = XVXT + Xtemp.transpose() * Xtemp / sigma; //4x4
    YZBX = YZBX + Xtemp.transpose() * (Ytemp - Ztemp * bVec) / sigma; // 4x1
    
  }
  
  
  
  beta1New = XVXT.inverse() * YZBX;
  
  XVXT = Eigen::MatrixXd::Zero(p1, p1);
  YZBX = Eigen::MatrixXd::Zero(p1, 1);
  
  for (int i = 0; i < numSubj; i++) {
    
    Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
    Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
    Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
    Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
    
    
    Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[1]);
    Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[1]);
    Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[1]);
    
    double sigma = Rcpp::as<double>(sigmaList[1]);
    Eigen::MatrixXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
    
    XVXT = XVXT + Xtemp.transpose() * Xtemp / sigma; //4x4
    YZBX = YZBX + Xtemp.transpose() * (Ytemp - Ztemp * bVec) / sigma; // 4x1
    
    
    
  }
  
  beta2New = XVXT.inverse() * YZBX;
  
  betaNew << beta1New, beta2New;
  Eigen::MatrixXd FUNB = Eigen::MatrixXd::Zero(8, numSubj);
  
  /*******************
   *  SIG CALCULATION
   ********************/
  int p1a = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(ZList[0])[0]).cols();
  int p2a = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(ZList[0])[1]).cols();
  
  Eigen::MatrixXd SigNew = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
  Eigen::MatrixXd num = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
  Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(p1a + p2a);
  
  int count = 0;
  for (int i = 0; i < numSubj; i++) {
    Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
    Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    
    Eigen::VectorXd bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
    Eigen::VectorXd bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
    
    bVeci << bVec1, bVec2;
    
    num += sigmai + MultVVoutprod(bVeci);
    
  }
  
  SigNew = num / numSubj;
  
  
  /*******************
   *  SIGMA CALCULATION
   ********************/
  
  // quadriture
  int nijSum = 0;
  double mu = 0;
  double tau;
  double tausq;
  double num1q = 0;
  double num2q = 0;
  double sigma1q;
  double sigma2q;
  
  for (int i = 0; i < numSubj; i++) {
    Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
    Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
    Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
    Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
    Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[0]);
    Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[0]);
    Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[0]);
    Eigen::MatrixXd bVeci = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
    Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    Eigen::MatrixXd sigigg = sigmai.block(0, 0, p1a, p1a);
    
    Rcpp::List mdataList = Rcpp::as<Rcpp::List>(mdata[0]);
    int numRep = Rcpp::as<int>(mdataList[i]);
    
    for (int nij = 0; nij < numRep; nij++) {
      
      mu = MultVV(Xtemp.row(nij), beta1New) + MultVV(Ztemp.row(nij), bVeci);
      tausq = Ztemp.row(nij) * sigigg * Ztemp.row(nij).transpose();
      tau = sqrt(tausq);
      
      //std::cout << "Mu1: " << mu << std::endl;
      //std::cout << "Tau-squared1: " << tausq << std::endl;
      //std::cout << "Tau1: " << tau << std::endl;
      
      for (int c = 0; c < numWeight; c++) {
        num1q += weight(c) * pow((Ytemp(nij, 0) - mu - tau * absc(c)), 2);
      }
    }
    
    nijSum += numRep;
    
  }
  
  // std::cout << "Tau-squared1: " << tausq << std::endl;
  // std::cout << "Mu1: " << mu << std::endl;
  // std::cout << "num1 " << num1q << std::endl;
  
  for (int i = 0; i < numSubj; i++) {
    Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
    Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
    Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
    Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
    Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[1]);
    Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[1]);
    Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[1]);
    Eigen::MatrixXd bVeci = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
    Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    Eigen::MatrixXd sigigg = sigmai.block(p1a, p1a, p2a, p2a);
    
    Rcpp::List mdataList = Rcpp::as<Rcpp::List>(mdata[1]);
    int numRep = Rcpp::as<int>(mdataList[i]);
    
    for (int nij = 0; nij < numRep; nij++) {
      mu = MultVV(Xtemp.row(nij), beta2New) + MultVV(Ztemp.row(nij), bVeci);
      
      tausq = Ztemp.row(nij) * sigigg * Ztemp.row(nij).transpose();//MultVVoutprod(Ztemp.row(nij))*sigigg;
      tau = sqrt(tausq);
      
      //std::cout << "Mu2: " << mu << std::endl;
      //std::cout << "Tau-squared2: " << tausq << std::endl;
      //
      
      for (int c = 0; c < numWeight; c++) {
        num2q += weight(c) * pow((Ytemp(nij, 0) - mu - tau * absc(c)), 2);
      }
    }
    
  }
  
  // std::cout << "Tau-squared2: " << tausq << std::endl;
  // std::cout << "Mu2: " << mu << std::endl;
  // std::cout << "num2 " << num2q << std::endl;
  // 
  // std::cout << "nijSum " << nijSum << std::endl;
  
  sigma1q = num1q / nijSum;
  sigma2q = num2q / nijSum;
  
  
  /*******************
   *  HAZ CALCULATION
   ********************/
  
  double denomH = 0;
  
  double dem1q = 0;
  double dem2q = 0;
  int a = H01.rows();
  int b = H02.rows();
  int risk1_index = a - 1;
  int risk2_index = b - 1;
  
  for (int i = 0; i < numSubj; i++) {
    Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
    Rcpp::List alphaListElement1 = Rcpp::as<Rcpp::List>(alphaList[0]); // get first risk
    Rcpp::List alphaListElement2 = Rcpp::as<Rcpp::List>(alphaList[1]); // get second risk
    Eigen::VectorXd alpha11 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[0]); // get alpha1// risk 1b1
    Eigen::VectorXd alpha12 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[1]); // get alpha2// risk 1b2
    Eigen::VectorXd alpha21 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[0]); // get alpha1// risk 2b1
    Eigen::VectorXd alpha22 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[1]); // get alpha2// risk 2b2
    Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    
    Eigen::VectorXd alpha1 = Eigen::VectorXd::Zero(p1a + p2a);
    Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(p1a + p2a);
    Eigen::VectorXd bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
    Eigen::VectorXd bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
    
    alpha1 << alpha11, alpha12; // diff bio, same alpha
    bVeci << bVec1, bVec2;
    double muH1, tausq1, tau1;
    muH1 = MultVV(W.row(i), gamma1) + MultVV(alpha1, bVeci);
    tausq1 = alpha1.transpose() * sigmai * alpha1;
    tau1 = sqrt(tausq1);
    
    for (int c = 0; c < numWeight; c++) {
      dem1q += weight(c) * exp(muH1 + tau1 * absc(c));
    }
    
    if (cmprsk(i) == 1) {
      
      // check last subject
      if (i == numSubj - 1)
      {
        H01(risk1_index, 2) = H01(risk1_index, 1) / dem1q;
        risk1_index--;
      }
      // check if time change
      else if (survtime(i + 1) != survtime(i))
      {
        H01(risk1_index, 2) = H01(risk1_index, 1) / dem1q;
        risk1_index--;
      }
      // every other subject
      else
      {
        for (i = i + 1; i < numSubj; i++)
        {
          
          bListElement = Rcpp::as<Rcpp::List>(bList[i]);
          bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
          bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
          bVeci << bVec1, bVec2;
          sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
          muH1 = MultVV(W.row(i), gamma1) + MultVV(alpha1, bVeci);
          tausq1 = alpha1.transpose() * sigmai * alpha1;
          tau1 = sqrt(tausq1);
          for (int c = 0; c < numWeight; c++) {
            dem1q += weight(c) * exp(muH1 + tau1 * absc(c));
          }
          
          if (i == numSubj - 1)
          {
            H01(risk1_index, 2) = H01(risk1_index, 1) / dem1q;
            risk1_index--;
            break;
          }
          else if (survtime(i + 1) != survtime(i))
          {
            H01(risk1_index, 2) = H01(risk1_index, 1) / dem1q;
            risk1_index--;
            break;
          }
          else continue;
        }
      }
      
    }
    else continue;
    
    
  }
  
  for (int i = 0; i < numSubj; i++) {
    Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
    Rcpp::List alphaListElement1 = Rcpp::as<Rcpp::List>(alphaList[0]); // get first risk
    Rcpp::List alphaListElement2 = Rcpp::as<Rcpp::List>(alphaList[1]); // get second risk
    Eigen::VectorXd alpha11 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[0]); // get alpha1// ris k 1 b1
    Eigen::VectorXd alpha12 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[1]); // get alpha2// risk 1b2
    Eigen::VectorXd alpha21 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[0]); // get alpha1// risk 2b1
    Eigen::VectorXd alpha22 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[1]); // get alpha2// risk 2b2
    Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    
    Eigen::VectorXd alpha2 = Eigen::VectorXd::Zero(p1a + p2a);
    Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(p1a + p2a);
    Eigen::VectorXd bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
    Eigen::VectorXd bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
    
    alpha2 << alpha21, alpha22; //diff bio, same alpha
    bVeci << bVec1, bVec2;
    double muH2, tausq2, tau2;
    muH2 = MultVV(W.row(i), gamma2) + MultVV(alpha2, bVeci);
    tausq2 = alpha2.transpose() * sigmai * alpha2;
    tau2 = sqrt(tausq2);
    
    for (int c = 0; c < numWeight; c++) {
      dem2q += weight(c) * exp(muH2 + tau2 * absc(c));
    }
    
    
    if (cmprsk(i) == 2) {
      
      // check last subject
      if (i == numSubj - 1)
      {
        H02(risk2_index, 2) = H02(risk2_index, 1) / dem2q;
        risk2_index--;
      }
      // check if time change
      else if (survtime(i + 1) != survtime(i))
      {
        H02(risk2_index, 2) = H02(risk2_index, 1) / dem2q;
        risk2_index--;
      }
      // every other subject
      else
      {
        for (i = i + 1; i < numSubj; i++)
        {
          bListElement = Rcpp::as<Rcpp::List>(bList[i]);
          bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
          bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
          bVeci << bVec1, bVec2;
          muH2 = MultVV(W.row(i), gamma2) + MultVV(alpha2, bVeci);
          sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
          tausq2 = alpha2.transpose() * sigmai * alpha2;
          tau2 = sqrt(tausq2);
          
          for (int c = 0; c < numWeight; c++) {
            dem2q += weight(c) * exp(muH2 + tau2 * absc(c));
          }
          
          if (i == numSubj - 1)
          {
            H02(risk2_index, 2) = H02(risk2_index, 1) / dem2q;
            risk2_index--;
            break;
          }
          else if (survtime(i + 1) != survtime(i))
          {
            H02(risk2_index, 2) = H02(risk2_index, 1) / dem2q;
            risk2_index--;
            break;
          }
          else continue;
        }
      }
      
    }
    else continue;
    
  }
  
  /*******************
   *  PHI calculation - Quadriture
   ********************/	
  
  double scalef;
  double scalefH01 = 0;
  double scalefH02 = 0;
  
  risk1_index = a-1;
  risk2_index = b-1;
  
  int dimW = gamma1.size();
  
  Eigen::VectorXd scalel = Eigen::VectorXd::Zero(p1a + p2a);
  Eigen::MatrixXd scalellT = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
  Eigen::MatrixXd scalewl = Eigen::MatrixXd::Zero(dimW, p1a + p2a);
  
  Eigen::VectorXd Sw_new = Eigen::VectorXd::Zero(dimW);
  Eigen::VectorXd Sw_inter = Eigen::VectorXd::Zero(dimW);
  Eigen::MatrixXd Sww_new = Eigen::MatrixXd::Zero(dimW, dimW);
  Eigen::VectorXd Sl_new = Eigen::VectorXd::Zero(p1a + p2a);
  Eigen::VectorXd Sl_inter = Eigen::VectorXd::Zero(p1a + p2a);
  Eigen::MatrixXd Sll_new = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
  Eigen::MatrixXd Swl_new = Eigen::MatrixXd::Zero(dimW, p1a + p2a);
  
  Eigen::MatrixXd  wwT = Eigen::MatrixXd::Zero(dimW, dimW);
  Eigen::MatrixXd  llT = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
  Eigen::MatrixXd  SwwT = Eigen::MatrixXd::Zero(dimW, dimW);
  Eigen::MatrixXd  SllT = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
  Eigen::VectorXd  Sw = Eigen::VectorXd::Zero(dimW);
  Eigen::VectorXd  Sl = Eigen::VectorXd::Zero(p1a + p2a);
  Eigen::MatrixXd Swl = Eigen::MatrixXd::Zero(dimW, p1a + p2a);
  
  Rcpp::List alphaListElement1 = Rcpp::as<Rcpp::List>(alphaList[0]); // get first risk
  Rcpp::List alphaListElement2 = Rcpp::as<Rcpp::List>(alphaList[1]); // get second risk
  Eigen::VectorXd alpha11 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[0]); // get alpha1// risk 1 b1
  Eigen::VectorXd alpha12 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[1]); // get alpha2// risk 1 b2
  Eigen::VectorXd alpha21 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[0]); // get alpha1// risk 2 b1
  Eigen::VectorXd alpha22 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[1]); // get alpha2// risk 2 b2
  
  Eigen::VectorXd alpha1 = Eigen::VectorXd::Zero(p1a + p2a);
  Eigen::VectorXd alpha2 = Eigen::VectorXd::Zero(p1a + p2a);
  alpha1 << alpha11, alpha12; // diff bio, same alpha
  alpha2 << alpha21, alpha22; // diff bio, same alpha
  
  
  for (int i = 0; i < numSubj; i++) {
    
    scalef = 0;
    scalel = Eigen::VectorXd::Zero(p1a + p2a);
    scalellT = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
    scalewl = Eigen::MatrixXd::Zero(dimW, p1a + p2a);
    
    Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
    
    Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(p1a + p2a);
    Eigen::VectorXd bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
    Eigen::VectorXd bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
    
    bVeci << bVec1, bVec2;
    
    
    // Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    Eigen::MatrixXd BAssociation = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    Eigen::VectorXd latent = bVeci;
    double mu1, tausq, tau;
    tausq = alpha1.transpose() * BAssociation * alpha1;
    tau = sqrt(tausq);
    
    Eigen::VectorXd w = W.row(i);
    mu1 = MultVV(w, gamma1) + MultVV(alpha1, bVeci);
    wwT = MultVVoutprod(w);
    
    // quadriture part/scale
    for (int c = 0; c < numWeight; c++) {
      scalef += weight(c) * exp(mu1 + tau * absc(c));
      scalel += weight(c) * exp(mu1 + tau * absc(c)) * (latent + BAssociation * alpha1 * absc(c) / tau);
      scalellT += weight(c) * exp(mu1 + tau * absc(c)) * (MultVVoutprod(latent + BAssociation*alpha1*absc(c)/tau) + (-MultVVoutprod(BAssociation*alpha1)/pow(tau,3)+BAssociation/tau)*absc(c));
      scalewl += weight(c) * exp(mu1 + tau * absc(c)) * w * (latent + BAssociation*alpha1*absc(c)/tau).transpose();
    }
    
    
    // for S
    w *= scalef; //sum_c w_c exp(mu + tau*v) * W
    Sw += w; //sum_i sum_c w_c exp(mu + tau*v) * W
    Sl += scalel; // sum_i w_c * exp(mu+tau*v) * (M+B*alpha*v/tau)
    
    // for I
    wwT *= scalef; //exp(mu+tau)wwT
    SwwT += wwT; //I_gamma //sum exp(mu)bbT
    SllT += scalellT; //I_alpha // sum_i sum_c w_c*exp(mu+tau*v) * ((M+B*alpha*v/tau)(M+B*alpha*v/tau)T+(-(B*alpha1)(B*alpha1)^T/tau^3+B/tau)*v)
    Swl += scalewl;// I_ga // sum_i sum_c w_c * exp(mu + tau*v) * w * (M+B*alpha*v/tau)^T
    
    if (cmprsk(i) == 1)
    {
      
      scalefH01 = H01(risk1_index, 2);
      
      // i
      Sww_new += scalefH01*SwwT;
      Sll_new += scalefH01*SllT;
      Swl_new += scalefH01*Swl;
      
      // s
      Sw_new += scalefH01*Sw;
      Sl_new += scalefH01*Sl;
      
      risk1_index--;
      
    }
    
    
  }
  
  
  
  for (int i = 0; i < numSubj; i++)
  {
    Eigen::VectorXd alpha1 = Eigen::VectorXd::Zero(p1a + p2a);
    Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(p1a + p2a);
    Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
    Eigen::VectorXd bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
    Eigen::VectorXd bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
    bVeci << bVec1, bVec2;
    Eigen::VectorXd latent = bVeci;
    
    if (cmprsk(i) == 1) {
      Sw_inter += W.row(i);
      Sl_inter += latent;
    }
  }
  
  //NR update
  
  Eigen::VectorXd Sfull_inter = Eigen::VectorXd::Zero(dimW + p1a + p2a);
  Eigen::VectorXd Sfull_new = Eigen::VectorXd::Zero(dimW + p1a + p2a);
  Eigen::MatrixXd info = Eigen::MatrixXd::Zero(dimW + p1a + p2a, dimW + p1a + p2a);
  
  Sfull_inter << Sw_inter, Sl_inter;
  Sfull_new << Sw_new, Sl_new;
  
  
  // start row, start column, how many rows, how many col
  info.block(0, 0, dimW, dimW) = Sww_new;
  info.block(0, dimW, dimW, p1a + p2a) = Swl_new;
  info.block(dimW, 0, p1a + p2a, dimW) = Swl_new.transpose();
  info.block(dimW, dimW, p1a + p2a, p1a + p2a) = Sll_new;
  
  // NR update
  Eigen::VectorXd phi1q = Eigen::VectorXd::Zero(dimW + p1a + p2a);
  phi1q << gamma1, alpha1;
  phi1q += info.inverse() * (Sfull_inter - Sfull_new);
  
  /////////
  // PHI2
  /////////
  
  dimW = gamma2.size();
  Sw_new = Eigen::VectorXd::Zero(dimW);
  Sw_inter = Eigen::VectorXd::Zero(dimW);
  Sww_new = Eigen::MatrixXd::Zero(dimW, dimW);
  Sl_new = Eigen::VectorXd::Zero(p1a + p2a);
  Sl_inter = Eigen::VectorXd::Zero(p1a + p2a);
  Sll_new = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
  Swl_new = Eigen::MatrixXd::Zero(dimW, p1a + p2a);
  
  wwT = Eigen::MatrixXd::Zero(dimW, dimW);
  llT = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
  SwwT = Eigen::MatrixXd::Zero(dimW, dimW);
  SllT = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
  Sw = Eigen::VectorXd::Zero(dimW);
  Sl = Eigen::VectorXd::Zero(p1a + p2a);
  Swl = Eigen::MatrixXd::Zero(dimW, p1a + p2a);
  
  
  for (int i = 0; i < numSubj; i++) {
    
    scalef = 0;
    scalel = Eigen::VectorXd::Zero(p1a + p2a);
    scalellT = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
    scalewl = Eigen::MatrixXd::Zero(dimW, p1a + p2a);
    
    Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
    
    Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(p1a + p2a);
    Eigen::VectorXd bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
    Eigen::VectorXd bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
    
    bVeci << bVec1, bVec2;
    
    
    // Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    Eigen::MatrixXd BAssociation = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
    Eigen::VectorXd latent = bVeci;
    double mu2, tausq, tau;
    tausq = alpha2.transpose() * BAssociation * alpha2;
    tau = sqrt(tausq);
    
    Eigen::VectorXd w = W.row(i);
    mu2 = MultVV(w, gamma2) + MultVV(alpha2, latent);
    wwT = MultVVoutprod(w);
    
    for (int c = 0; c < numWeight; c++) {
      scalef += weight(c) * exp(mu2 + tau * absc(c)); // Sgamma
      scalel += weight(c) * exp(mu2 + tau * absc(c)) * (latent + BAssociation * alpha2 * absc(c) / tau); //S_alpha
      scalellT += weight(c) * exp(mu2 + tau * absc(c)) * (MultVVoutprod(latent + BAssociation*alpha2*absc(c)/tau) + (-MultVVoutprod(BAssociation*alpha2)/pow(tau,3)+BAssociation/tau)*absc(c)); //I_alpha
      scalewl += weight(c) * exp(mu2 + tau * absc(c)) * w * (latent + BAssociation*alpha2*absc(c)/tau).transpose(); //I_ga
    }
    
    // for S
    w *= scalef; //sum_c w_c exp(mu + tau*v) * W
    Sw += w; //sum_i sum_c w_c exp(mu + tau*v) * W
    Sl += scalel; // sum_i w_c * exp(mu+tau*v) * (M+B*alpha*v/tau)
    
    // for I
    wwT *= scalef; //exp(mu+tau)wwT
    SwwT += wwT; //I_gamma //sum exp(mu)bbT
    SllT += scalellT; //I_alpha // sum_i sum_c w_c*exp(mu+tau*v) * ((M+B*alpha*v/tau)(M+B*alpha*v/tau)T+(-(B*alpha1)(B*alpha1)^T/tau^3+B/tau)*v)
    Swl += scalewl;// I_ga // sum_i sum_c w_c * exp(mu + tau*v) * w * (M+B*alpha*v/tau)^T
    
    
    if (cmprsk(i) == 2)
    {
      
      scalefH02 = H02(risk2_index, 2);
      
      // i
      Sww_new += scalefH02*SwwT;
      Sll_new += scalefH02*SllT;
      Swl_new += scalefH02*Swl;
      
      // s
      Sw_new += scalefH02*Sw;
      Sl_new += scalefH02*Sl;
      
      risk2_index--;
      
    }
    
  }
  
  
  for (int i = 0; i < numSubj; i++)
  {
    Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(p1a + p2a);
    Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
    Eigen::VectorXd bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
    Eigen::VectorXd bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
    bVeci << bVec1, bVec2;
    Eigen::VectorXd latent = bVeci;
    
    if (cmprsk(i) == 2) {
      Sw_inter += W.row(i);
      Sl_inter += latent;
    }
  }
  
  
  //NR update
  
  Sfull_inter = Eigen::VectorXd::Zero(dimW + p1a + p2a);
  Sfull_new = Eigen::VectorXd::Zero(dimW + p1a + p2a);
  info = Eigen::MatrixXd::Zero(dimW + p1a + p2a, dimW + p1a + p2a);
  
  Sfull_inter << Sw_inter, Sl_inter;
  Sfull_new << Sw_new, Sl_new;
  
  // start row, start column, how many rows, how many col
  info.block(0, 0, dimW, dimW) = Sww_new;
  info.block(0, dimW, dimW, p1a + p2a) = Swl_new;
  info.block(dimW, 0, p1a + p2a, dimW) = Swl_new.transpose();
  info.block(dimW, dimW, p1a + p2a, p1a + p2a) = Sll_new;
  
  // NR update
  Eigen::VectorXd phi2q = Eigen::VectorXd::Zero(dimW + p1a + p2a);
  phi2q << gamma2, alpha2;
  phi2q += info.inverse() * (Sfull_inter - Sfull_new); //+I * S
  
  
  return Rcpp::List::create(Rcpp::Named("FUNB") = FUNB, Rcpp::Named("beta1") = beta1New, Rcpp::Named("beta2") = beta2New, 
                                        Rcpp::Named("beta") = betaNew, Rcpp::Named("Sig") = SigNew, Rcpp::Named("sigma1") = sigma1q, Rcpp::Named("sigma2") = sigma2q, 
                                                    Rcpp::Named("H01") = H01, Rcpp::Named("H02") = H02, 
                                                    Rcpp::Named("phi1") = phi1q, Rcpp::Named("phi2") = phi2q);
  
}


