#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]


Rcpp::List getMNA(Rcpp::List XList, Rcpp::List YList, Rcpp::List ZList, Eigen::MatrixXd& W,
	Rcpp::List mdata, Rcpp::List mdataSList, 
	Rcpp::List bList, Rcpp::List sigmaList, Rcpp::List sigmaiList,
	Eigen::VectorXd weight, Eigen::VectorXd absc,
	Eigen::MatrixXd H01, Eigen::MatrixXd H02, Eigen::VectorXd& survtime, Eigen::VectorXd cmprsk,
	Eigen::VectorXd& gamma1, Eigen::VectorXd& gamma2, Rcpp::List alphaList, const Eigen::MatrixXd& xsmatrix,
	const Eigen::MatrixXd& wsmatrix,
	const Eigen::VectorXd& CUH01,
	const Eigen::VectorXd& CUH02,
	const Eigen::VectorXd& HAZ01,
	const Eigen::VectorXd& HAZ02, const Eigen::MatrixXd& Sig,
	Rcpp::List betaList){


	Eigen::MatrixXd H01q = H01;
	Eigen::MatrixXd H02q = H02;
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



for (int i = 0; i < numSubj; i++) {
			// Ensure that the current subject exists in XList, YList, ZList
			//if (i < XList.size() && i < YList.size() && i < ZList.size()) {
			Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
			Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
			Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
			Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
	
			//for (int g = 0; g < numBio; g++) {
	
			//        Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[g]);
			//        Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[g]);
			//        Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[g]);
	
			//        double sigma = Rcpp::as<double>(sigmaList[0]);
			//        // Convert bList[i][g] to Eigen::VectorXd if it is not already
			//        Eigen::MatrixXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
	
			//        XVXT = XVXT + Xtemp.transpose() * Xtemp / sigma; //4x4
			//        YZBX = YZBX + Xtemp.transpose() * (Ytemp - Ztemp * bVec) / sigma; // 4x1
	
			//        Eigen::VectorXd beta = XVXT.inverse() * YZBX;
			//        betas.push_back(beta);
	
	
			//}
	
			Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[0]);
			Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[0]);
			Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[0]);
	
			double sigma = Rcpp::as<double>(sigmaList[0]);
	
	
			//std::cout << "X\n" << Xtemp << std::endl;
			//std::cout << "Y\n" << Ytemp << std::endl;
			//std::cout << "Z\n" << Ztemp << std::endl;
	
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

	return Rcpp::List::create(Rcpp::Named("FUNB") = FUNB, Rcpp::Named("beta1") = beta1New, Rcpp::Named("beta2") = beta2New, 
		Rcpp::Named("beta") = betaNew);
		//	Rcpp::Named("beta1") = beta1New, Rcpp::Named("beta2") = beta2New, Rcpp::Named("beta") = betaNew,
		//	Rcpp::Named("Sig") = SigNew, //Rcpp::Named("b") = bVeci,
		//	Rcpp::Named("sigma1") = sigma1, Rcpp::Named("sigma2") = sigma2, Rcpp::Named("sigma1q") = sigma1q, Rcpp::Named("sigma2q") = sigma2q,
		//	Rcpp::Named("H01") = H01, Rcpp::Named("H02") = H02, Rcpp::Named("H01q") = H01q, Rcpp::Named("H02q") = H02q,
		//		Rcpp::Named("phi1") = phi1, Rcpp::Named("phi2") = phi2);// , //Rcpp::Named("phi2q") = phi2q,
				//Rcpp::Named("check") = check, Rcpp::Named("check2") = check2,
				//Rcpp::Named("check3") = check3, Rcpp::Named("check4") = check4);
//return Rcpp::List::create(
//	Rcpp::Named("beta1") = beta1New, Rcpp::Named("beta2") = beta2New, Rcpp::Named("beta") = betaNew,
//	Rcpp::Named("Sig") = SigNew, //Rcpp::Named("b") = bVeci,
//	Rcpp::Named("sigma1") = sigma1, Rcpp::Named("sigma2") = sigma2, Rcpp::Named("sigma1q") = sigma1q, Rcpp::Named("sigma2q") = sigma2q,
//	Rcpp::Named("H01") = H01, Rcpp::Named("H02") = H02, Rcpp::Named("H01q") = H01q, Rcpp::Named("H02q") = H02q,
//		Rcpp::Named("phi1") = phi1, Rcpp::Named("phi2") = phi2);// , //Rcpp::Named("phi2q") = phi2q,
		//Rcpp::Named("check") = check, Rcpp::Named("check2") = check2,
		//Rcpp::Named("check3") = check3, Rcpp::Named("check4") = check4);

}


