#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]


Rcpp::List getMaGH(Rcpp::List XList, Rcpp::List YList, Rcpp::List ZList, Eigen::MatrixXd& W,
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

	/*******************
*  ADAPT-GAUSS BETA CALCULATION
********************/
	Rcpp::List mdataList = Rcpp::as<Rcpp::List>(mdata[0]);
	//int numSubj = mdataList.size();
	int p1a = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(ZList[0])[0]).cols();
	int p2a = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(ZList[0])[1]).cols();

	double dem, cuh01, cuh02, haz01, haz02, xgamma1, xgamma2, temp, mu, zb;
	Eigen::VectorXd c(p1a + p2a);
	Eigen::VectorXd bi(p1a + p2a);
	Eigen::VectorXd bi2(p1a + p2a);
	Eigen::VectorXd weightbi(p1a + p2a);
	Eigen::VectorXd ri(p1a + p2a);
	Eigen::VectorXd rii(p1a + p2a);
	Eigen::MatrixXd Hi(p1a + p2a, p1a + p2a);
	Eigen::MatrixXd Hi2(p1a + p2a, p1a + p2a);

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig.inverse(), Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::VectorXd eigenSQ = svd.singularValues();
	//int i, j, q, t, db;
	for (int i = 0; i < eigenSQ.size(); i++) {
		eigenSQ(i) = sqrt(eigenSQ(i));
	}
	Eigen::MatrixXd SigSQRT = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();

	//std::cout << "Sigsqrtr\n " << SigSQRT << std::endl;
	Eigen::MatrixXd FUNB = Eigen::MatrixXd::Zero(p1a+p2a, numSubj);

	int point = wsmatrix.rows();
	int numRep;
	for (int i = 0; i < numSubj; i++)
	{
		dem = 0;
		numRep = mdataList(i);
		cuh01 = CUH01(i);
		cuh02 = CUH02(i);
		haz01 = HAZ01(i);
		haz02 = HAZ02(i);

		xgamma1 = MultVV(W.row(i), gamma1);
		xgamma2 = MultVV(W.row(i), gamma2);

		Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);

		Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(p1a + p2a);

		//Eigen::MatrixXd bVeci = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
		Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
		//Eigen::MatrixXd sigig = sigmai.block(0, 0, p1a, p1a);

		Eigen::VectorXd bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
		Eigen::VectorXd bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);

	
		bVeci << bVec1, bVec2;

		for (int ct = 0; ct < p1a+p2a; ct++) Hi.row(ct) = sigmai.row(ct);
        //Hi.row(ct) = Poscov.row(j * p1a + ct);
			
		//std::cout << "Hi\n " << Hi << std::endl;


		Eigen::JacobiSVD<Eigen::MatrixXd> svd(Hi, Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::VectorXd eigenSQ = svd.singularValues();
		for (int ct = 0; ct < eigenSQ.size(); ct++) {
			eigenSQ(ct) = sqrt(eigenSQ(ct));
		}

		Hi2 = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();

		//std::cout << "Hi2\n " << Hi2 << std::endl;

		for (int db = 0; db < point; db++) {
			c = xsmatrix.row(db);
			weightbi = wsmatrix.row(db);
			////bii = Posbi.row(j);
			 //4x1 4x4
			ri = bVeci + sqrt(2) * Hi2 * c;
			rii = SigSQRT * ri;
			//temp = exp(705);
            temp = 0;
            
			//biomarker; take first 2 values of ri. line 140-end of loop;; take first 2 ri, bottom 2 ri for each biomarker
			// fopr survval/prior dist of re do whole thing

			for (int ct = 0; ct < p1a+p2a; ct++) c(ct) = ri(ct);

			//std::cout << "c\n" << c << std::endl;
			//std::cout << "c\n" << ri << std::endl;

			Rcpp::List mdataList = Rcpp::as<Rcpp::List>(mdata[0]);
			int numRep = Rcpp::as<int>(mdataList[i]);
			Rcpp::List mdataS = Rcpp::as<Rcpp::List>(mdataSList[0]);

			

			for (int g = 0; g < numBio; g++) {

				Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
				Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
				Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
				//Rcpp::List betag = Rcpp::as<Rcpp::List>(BetaList[g]);
				Eigen::VectorXd betag = Rcpp::as<Eigen::VectorXd>(betaList[g]);
				Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[g]);
				Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[g]);
				Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[g]);

				double sigma = Rcpp::as<double>(sigmaList[g]);

				int index = 0;
				Eigen::VectorXd cbio = c.segment(index, Ztemp.cols());
				index += Ztemp.cols();
				
				for (int nij = 0; nij < numRep; nij++) {
					mu = MultVV(Xtemp.row(nij), betag);
					//std::cout << "xtemp\n" << Xtemp.row(nij) << std::endl;
					//std::cout << "beta" << beta1 << std::endl;

					//std::cout << "ind" << index << std::endl;
					//std::cout << "b" << bi << std::endl;

					zb = MultVV(Ztemp.row(nij), cbio);
					//double tempjr = temp;
					//temp *= exp(-1 / (2 * sigma) * pow((Ytemp(nij) - mu - zb), 2));
                    temp += -1 / (2 * sigma) * pow((Ytemp(nij) - mu - zb), 2);
            

				}

			}

			Rcpp::List alphaListElement1 = Rcpp::as<Rcpp::List>(alphaList[0]); // get first  risk
			Rcpp::List alphaListElement2 = Rcpp::as<Rcpp::List>(alphaList[1]); // get second risk
			Eigen::VectorXd alpha11 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[0]); // get alpha1// ris k 1 b1
			Eigen::VectorXd alpha12 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[1]); // get alpha2// risk 1b2
			Eigen::VectorXd alpha21 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[0]); // get alpha1// risk 2b1
			Eigen::VectorXd alpha22 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[1]); // get alpha2// risk 2b2
			//Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);

			Eigen::VectorXd alpha1 = Eigen::VectorXd::Zero(p1a + p2a);
			Eigen::VectorXd alpha2 = Eigen::VectorXd::Zero(p1a + p2a);

			alpha1 << alpha11, alpha12; // diff bio, same alpha
			alpha2 << alpha21, alpha22; // diff bio, same alpha
			

			//if (cmprsk(i) == 1)  temp *= haz01 * exp(xgamma1 + MultVV(alpha1, bi));
			//if (cmprsk(i) == 2)  temp *= haz02 * exp(xgamma2 + MultVV(alpha2, bi));
			if (cmprsk(i) == 1) { 
				temp += log(haz01) + (xgamma1 + MultVV(alpha1, c)); 
				//std::cout << "log" << log(haz01) + (xgamma1 + MultVV(alpha1, bi)) << std::endl;
			}
			if (cmprsk(i) == 2) {
				temp += log(haz02) + (xgamma2 + MultVV(alpha2, c));
				//std::cout << "log" << log(haz02) + (xgamma2 + MultVV(alpha2, bi)) << std::endl;
			}
			
			//std::cout << "temp2\n" << temp << std::endl;

			//temp *= exp(0 - cuh01 * exp(xgamma1 + MultVV(alpha1, bi)) - cuh02 * exp(xgamma2 + MultVV(alpha2, bi)));
            temp += - cuh01 * exp(xgamma1 + MultVV(alpha1, c)) - cuh02 * exp(xgamma2 + MultVV(alpha2, c));

			//std::cout << "temp3\n" << temp << std::endl;

			for (int ct = 0; ct < p1a+p2a; ct++) temp *= weightbi(ct);
			// prior distribution of RE; evaluvate whole vector
			bi2 = xsmatrix.row(db);
			//temp *= exp(-pow(rii.norm(), 2) / 2) * exp(pow(bi2.norm(), 2));
            temp += (-pow(rii.norm(), 2) / 2) + pow(bi2.norm(), 2);
                                                    
            temp = exp(temp);
            
			dem += temp;
			//calculate h(bi)
			FUNB.col(i) += temp * c;

			

		}                                                            
                                  
        FUNB.col(i)/=dem;
	}
	
	for (int i = 0; i < numSubj; i++) {
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

		YZBX = YZBX + Xtemp.transpose() * (Ytemp - Ztemp * FUNB.col(i)) / sigma; // 4x1

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
		YZBX = YZBX + Xtemp.transpose() * (Ytemp - Ztemp * FUNB.col(i)) / sigma; // 4x1



	}

	beta2New = XVXT.inverse() * YZBX;

	betaNew << beta1New, beta2New;

//for (int i = 0; i < numSubj; i++) {
//			// Ensure that the current subject exists in XList, YList, ZList
//			//if (i < XList.size() && i < YList.size() && i < ZList.size()) {
//			Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
//			Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
//			Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
//			Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
//	
//			//for (int g = 0; g < numBio; g++) {
//	
//			//        Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[g]);
//			//        Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[g]);
//			//        Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[g]);
//	
//			//        double sigma = Rcpp::as<double>(sigmaList[0]);
//			//        // Convert bList[i][g] to Eigen::VectorXd if it is not already
//			//        Eigen::MatrixXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
//	
//			//        XVXT = XVXT + Xtemp.transpose() * Xtemp / sigma; //4x4
//			//        YZBX = YZBX + Xtemp.transpose() * (Ytemp - Ztemp * bVec) / sigma; // 4x1
//	
//			//        Eigen::VectorXd beta = XVXT.inverse() * YZBX;
//			//        betas.push_back(beta);
//	
//	
//			//}
//	
//			Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[0]);
//			Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[0]);
//			Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[0]);
//	
//			double sigma = Rcpp::as<double>(sigmaList[0]);
//	
//	
//			//std::cout << "X\n" << Xtemp << std::endl;
//			//std::cout << "Y\n" << Ytemp << std::endl;
//			//std::cout << "Z\n" << Ztemp << std::endl;
//	
//			Eigen::MatrixXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
//	
//			XVXT = XVXT + Xtemp.transpose() * Xtemp / sigma; //4x4
//			YZBX = YZBX + Xtemp.transpose() * (Ytemp - Ztemp * bVec) / sigma; // 4x1
//
//		}
//	
//	
//	
//		beta1New = XVXT.inverse() * YZBX;
//	
//		XVXT = Eigen::MatrixXd::Zero(p1, p1);
//		YZBX = Eigen::MatrixXd::Zero(p1, 1);
//	
//		for (int i = 0; i < numSubj; i++) {
//	
//			Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
//			Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
//			Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
//			Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
//	
//	
//			Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[1]);
//			Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[1]);
//			Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[1]);
//	
//			double sigma = Rcpp::as<double>(sigmaList[1]);
//			Eigen::MatrixXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
//	
//			XVXT = XVXT + Xtemp.transpose() * Xtemp / sigma; //4x4
//			YZBX = YZBX + Xtemp.transpose() * (Ytemp - Ztemp * bVec) / sigma; // 4x1
//	
//			
//	
//		}
//	
//		beta2New = XVXT.inverse() * YZBX;
//	
//		betaNew << beta1New, beta2New;
//		Eigen::MatrixXd FUNB = Eigen::MatrixXd::Zero(8, numSubj);

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


