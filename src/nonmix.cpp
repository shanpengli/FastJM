#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]


Rcpp::List getNoQuad(Rcpp::List XList, Rcpp::List YList, Rcpp::List ZList, Eigen::MatrixXd& W,
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


		// ----------------------
		//   sigma
		// ----------------------

		double sigmai1 = 0;
		double num1 = 0;
		double sigmai2 = 0;
		double num2 = 0;
		double sigma1, sigma2;
		int p1a = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(ZList[0])[0]).cols();
		int p2a = Rcpp::as<Eigen::MatrixXd>(Rcpp::as<Rcpp::List>(ZList[0])[1]).cols();

		Eigen::MatrixXd ZZT = Eigen::MatrixXd::Zero(p1a, p1a);
		int nijSum = 0;

		// --- Normal Approximation ----------------------
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
			Eigen::MatrixXd sigig = sigmai.block(0, 0, p1a, p1a);
			Rcpp::List mdataList = Rcpp::as<Rcpp::List>(mdata[0]);
			int numRep = Rcpp::as<int>(mdataList[i]);

			for (int nij = 0; nij < numRep; nij++) {
				//double Ynij = Ytemp.row(nij)
				double epsilon = Ytemp(nij, 0) - MultVV(Xtemp.row(nij), beta1New);
	
				double zb = MultVV(Ztemp.row(nij), bVeci);
				Eigen::MatrixXd ZZT = MultVVoutprod(Ztemp.row(nij));
				Eigen::MatrixXd bbT = MultVVoutprod(bVeci);

				num1 += pow(epsilon, 2) - 2 * epsilon * zb + (ZZT * (sigig + bbT)).trace();
			}

			nijSum += numRep;
		}

	Eigen::MatrixXd SigE = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
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

	SigE = num / numSubj;

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
			int q = Rcpp::as<int>(mdataList[i]);

			for (int nij = 0; nij < q; nij++) {
				double epsilon = Ytemp(nij, 0) - MultVV(Xtemp.row(nij), beta2New);
				double zb = MultVV(Ztemp.row(nij), bVeci);
				Eigen::MatrixXd ZZT = MultVVoutprod(Ztemp.row(nij));
				Eigen::MatrixXd bbT = MultVVoutprod(bVeci);
				num2 += pow(epsilon, 2) - 2 * epsilon * zb + (ZZT * (sigigg + bbT)).trace();
			}

		}

		sigma1 = num1 / nijSum;
		sigma2 = num2 / nijSum;


		// ----------------------
		//   Hazard
		// ----------------------
		double denomH = 0;

		double dem1 = 0;
		double dem2 = 0;
		int a = H01.rows();
		int b = H02.rows();
		int risk1_index = a - 1;
		int risk2_index = b - 1;

		// // --- normal ----------------------
		for (int i = 0; i < numSubj; i++) {
			Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
			Rcpp::List alphaListElement1 = Rcpp::as<Rcpp::List>(alphaList[0]); // get first  risk
			Rcpp::List alphaListElement2 = Rcpp::as<Rcpp::List>(alphaList[1]); // get second risk
			Eigen::VectorXd alpha11 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[0]); // get alpha1// ris k 1 b1
			Eigen::VectorXd alpha12 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[1]); // get alpha2// risk 1b2
			Eigen::VectorXd alpha21 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[0]); // get alpha1// risk 2b1
			Eigen::VectorXd alpha22 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[1]); // get alpha2// risk 2b2
			Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);

			Eigen::VectorXd alpha1 = Eigen::VectorXd::Zero(p1a + p2a);
			Eigen::VectorXd alpha2 = Eigen::VectorXd::Zero(p1a + p2a);
			Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(p1a + p2a);
			Eigen::VectorXd bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
			Eigen::VectorXd bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);

			alpha1 << alpha11, alpha12; // diff bio, same alpha
			bVeci << bVec1, bVec2;
			double muH1, tau1, tausq1;
			muH1 = MultVV(W.row(i), gamma1) + MultVV(alpha1, bVeci);

			tausq1 = alpha1.transpose() * sigmai * alpha1;
			// tau1 = sqrt(tausq1);
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
						bListElement = Rcpp::as<Rcpp::List>(bList[i]);
						bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
						bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
						bVeci << bVec1, bVec2;
						sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
						muH1 = MultVV(W.row(i), gamma1) + MultVV(alpha1, bVeci);
						tausq1 = alpha1.transpose() * sigmai * alpha1;
						// tau1 = sqrt(tausq1);
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


		for (int i = 0; i < numSubj; i++) {
			Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
			Rcpp::List alphaListElement1 = Rcpp::as<Rcpp::List>(alphaList[0]); // get first risk
			Rcpp::List alphaListElement2 = Rcpp::as<Rcpp::List>(alphaList[1]); // get second risk
			Eigen::VectorXd alpha11 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[0]); // get alpha1// ris k 1 b1
			Eigen::VectorXd alpha12 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[1]); // get alpha2// risk 1b2
			Eigen::VectorXd alpha21 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[0]); // get alpha1// risk 2b1
			Eigen::VectorXd alpha22 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[1]); // get alpha2// risk 2b2
			Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);

			Eigen::VectorXd alpha1 = Eigen::VectorXd::Zero(p1a + p2a);
			Eigen::VectorXd alpha2 = Eigen::VectorXd::Zero(p1a + p2a);
			Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(p1a + p2a);
			Eigen::VectorXd bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
			Eigen::VectorXd bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);


			alpha2 << alpha21, alpha22; //diff bio, same alpha
			bVeci << bVec1, bVec2;

			double muH2, tausq2, tau2;
			muH2 = MultVV(W.row(i), gamma2) + MultVV(alpha2, bVeci);
			tausq2 = alpha2.transpose() * sigmai * alpha2;
			//tau2 = sqrt(tausq2);
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
						bListElement = Rcpp::as<Rcpp::List>(bList[i]);
						bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
						bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
						bVeci << bVec1, bVec2;
						sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
						muH2 = MultVV(W.row(i), gamma2) + MultVV(alpha2, bVeci);
						tausq2 = alpha2.transpose() * sigmai * alpha2;
						//tau2 = sqrt(tausq2);
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

		double scalefH01 = 0;
		double scalefH02 = 0;
		double scalef;

		risk1_index = a - 1;
		risk2_index = b - 1;

		int dimW = gamma1.size();

		Eigen::VectorXd Sw_new = Eigen::VectorXd::Zero(dimW);
		Eigen::VectorXd Sw_inter = Eigen::VectorXd::Zero(dimW);
		Eigen::MatrixXd Sww_new = Eigen::MatrixXd::Zero(dimW, dimW);
		Eigen::VectorXd Sl_new = Eigen::VectorXd::Zero(p1a + p2a);
		Eigen::VectorXd Sl_inter = Eigen::VectorXd::Zero(p1a + p2a);
		Eigen::MatrixXd Sll_new = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
		Eigen::MatrixXd Swl_new = Eigen::MatrixXd::Zero(dimW, p1a + p2a);

		Eigen::VectorXd latent;

		Eigen::MatrixXd  wwT = Eigen::MatrixXd::Zero(dimW, dimW);
		Eigen::MatrixXd  llT = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
		Eigen::MatrixXd  SwwT = Eigen::MatrixXd::Zero(dimW, dimW);
		Eigen::MatrixXd  SllT = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
		Eigen::VectorXd  Sw = Eigen::MatrixXd::Zero(dimW, dimW);
		Eigen::VectorXd  Sl = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
		Eigen::MatrixXd Swl = Eigen::MatrixXd::Zero(dimW, p1a + p2a);
		
		Rcpp::List alphaListElement1 = Rcpp::as<Rcpp::List>(alphaList[0]); // get first risk
		Rcpp::List alphaListElement2 = Rcpp::as<Rcpp::List>(alphaList[1]); // get second risk
		Eigen::VectorXd alpha11 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[0]); // get alpha1// risk 1 b1
		Eigen::VectorXd alpha12 = Rcpp::as<Eigen::VectorXd>(alphaListElement1[1]); // get alpha2// risk 1 b2
		Eigen::VectorXd alpha21 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[0]); // get alpha1// risk 2 b1
		Eigen::VectorXd alpha22 = Rcpp::as<Eigen::VectorXd>(alphaListElement2[1]); // get alpha2// risk 2 b2

		Eigen::VectorXd alpha1 = Eigen::VectorXd::Zero(p1a + p2a);
		alpha1 << alpha11, alpha12; // diff bio, same alpha


		for (int i = 0; i < numSubj; i++) {
			//scalef = 0;
			Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);

			Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(p1a + p2a);
			Eigen::VectorXd bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
			Eigen::VectorXd bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);


			bVeci << bVec1, bVec2;

			// Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
			Eigen::MatrixXd BAssociation = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
			Eigen::VectorXd latent = bVeci;
			double mu1, tausq;
			tausq = alpha1.transpose() * BAssociation * alpha1;

			Eigen::VectorXd w = W.row(i);
			//Eigen::VectorXd l = latent;
			Eigen::VectorXd l = BAssociation * alpha1 + bVeci;

			Eigen::MatrixXd wl = Eigen::MatrixXd::Zero(dimW, p1a + p2a);

			wl = w * l.transpose(); //iT

			mu1 = MultVV(w, gamma1) + MultVV(alpha1, bVeci);
			wwT = MultVVoutprod(w);
			//llT = MultVVoutprod(l);
			llT = MultVVoutprod(BAssociation * alpha1 + bVeci) + BAssociation;

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

				scalefH01 = H01(risk1_index, 2);
				//SXX *= scalefH01;
				// i
				SwwT *= scalefH01; //haz*exp(mu)wwT
				SllT *= scalefH01; //haz*mm * exp(mu)
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
				//     if (i == numSubj - 1)
				//     {
				//         scalefH01 = H01(risk1_index, 2);
				//         //SXX *= scalefH01;
				//         // i
				//         SwwT *= scalefH01; //haz*exp(mu)wwT
				//         SllT *= scalefH01; //haz*exp(mu) bbT
				//         Swl *= scalefH01;
				//         Sww_new += SwwT;
				//         Sll_new += SllT;
				//         Swl_new += Swl;

				//         SwwT /= scalefH01;
				//         SllT /= scalefH01;
				//         Swl /= scalefH01;

				//          // s
				//        Sw *= scalefH01;
				//        Sl *= scalefH01;
				//         Sw_new += Sw;
				//         Sl_new += Sl;
				//         Sw /= scalefH01;
				//         Sl /= scalefH01;


				//         risk1_index--;

				//     }

				//     else if (survtime(i + 1) != survtime(i))
				//     {
				//         scalefH01 = H01(risk1_index, 2);

				//         // i
				//         SwwT *= scalefH01;
				//         SllT *= scalefH01;
				//         Swl *= scalefH01;
				//         Sww_new += SwwT;
				//         Sll_new += SllT;
				//         Swl_new += Swl;
				//         SwwT /= scalefH01;
				//         SllT /= scalefH01;
				//         Swl /= scalefH01;

				//          // s
				//        Sw *= scalefH01;
				//        Sl *= scalefH01;
				//         Sw_new += Sw;
				//         Sl_new += Sl;
				//         Sw /= scalefH01;
				//         Sl /= scalefH01;
				//         risk1_index--;
				//     }
				//     else
				//     {
				//         for (i = i + 1; i < numSubj; i++)
				//         {

				//             bListElement = Rcpp::as<Rcpp::List>(bList[i]);
				//             bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
				//             bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);
				//             bVeci << bVec1, bVec2;
				//             latent = bVeci;
				//             w = W.row(i);
				//             l = latent;


				//             mu1 = MultVV(w, gamma1) + MultVV(alpha1, l);
				//             wwT =  MultVVoutprod(W.row(i));
				//             llT = MultVVoutprod(latent);

				//             scalef = exp(mu1);

				//             wl = Eigen::MatrixXd::Zero(dimW, p1a+p2a);

				//             wl = w * l.transpose();

				// // for I
				// wwT *= scalef;
				// llT *= scalef;
				// wl *= scalef;
				// SwwT += wwT;
				// SllT += llT;
				// Swl += wl;


				// // // for S
				// // wwT *= scalef;
				// // llT *= scalef;
				// // SwwT += wwT;
				// // SllT += llT;


				// // for S
				// w *= scalef;
				// l *= scalef;
				// Sw += w;
				// Sl += l;


				//             if (i == numSubj - 1)
				//             {
				//                 scalefH01 = H01(risk1_index, 2);
				//         SwwT *= scalefH01;
				//         SllT *= scalefH01;
				//         Swl *= scalefH01;
				//         Sww_new += SwwT;
				//         Sll_new += SllT;
				//         Swl_new += Swl;
				//         SwwT /= scalefH01;
				//         SllT /= scalefH01;
				//         Swl /= scalefH01;

				//          // s
				//        Sw *= scalefH01;
				//        Sl *= scalefH01;
				//         Sw_new += Sw;
				//         Sl_new += Sl;
				//         Sw /= scalefH01;
				//         Sl /= scalefH01;

				//         risk1_index--;
				//                 break;
				//             }
				//             else if (survtime(i + 1) != survtime(i))
				//             {
				//                 scalefH01 = H01(risk1_index, 2);

				//                 // i
				//         SwwT *= scalefH01;
				//         SllT *= scalefH01;
				//         Swl *= scalefH01;
				//         Sww_new += SwwT;
				//         Sll_new += SllT;
				//         Swl_new += Swl;
				//         SwwT /= scalefH01;
				//         SllT /= scalefH01;
				//         Swl /= scalefH01;

				//          // s
				//        Sw *= scalefH01;
				//        Sl *= scalefH01;
				//         Sw_new += Sw;
				//         Sl_new += Sl;
				//         Sw /= scalefH01;
				//         Sl /= scalefH01;
				//         risk1_index--;
				//                 break;
				//             }
				//             else continue;
				//         }
				//     }
				// }
				// else continue;
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
		Eigen::VectorXd phi1 = Eigen::VectorXd::Zero(dimW + p1a + p2a);
		phi1 << gamma1, alpha1;
		phi1 += info.inverse() * (Sfull_inter - Sfull_new);


		Sw_new = Eigen::VectorXd::Zero(dimW);
		Sw_inter = Eigen::VectorXd::Zero(dimW);
		Sww_new = Eigen::MatrixXd::Zero(dimW, dimW);
		Sl_new = Eigen::VectorXd::Zero(p1a + p2a);
		Sl_inter = Eigen::VectorXd::Zero(p1a + p2a);
		Sll_new = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
		Swl_new = Eigen::MatrixXd::Zero(dimW, p1a + p2a);
		// new matrix
		//Eigen::MatrixXd joint = Eigen::MatrixXd::Zero(p2 + p2, p2);
		// Eigen::MatrixXd  jjT = Eigen::MatrixXd::Zero(dimW   1a + p2a, dimW + p1a + p2a); // double check matrix dim might've switched
		// Eigen::MatrixXd SjjT = Eigen::MatrixXd::Zero(dimW + p1a + p2a, dimW + p1a + p2a);
		// Eigen::VectorXd Sj = Eigen::VectorXd::Zero(dimW + p1a + p2a);



		//     //MultVV(Z.row(j), FUNB.col(j) + tau * abscMat(0, j));

	//         //xb1 =  // fixed effect; contribute to LLH
	//         //zb = MultVV(FUNC, alpha1); // random effect component

		wwT = Eigen::MatrixXd::Zero(dimW, dimW);
		llT = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
		SwwT = Eigen::MatrixXd::Zero(dimW, dimW);
		SllT = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
		Sw = Eigen::MatrixXd::Zero(dimW, dimW);
		Sl = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
		Swl = Eigen::MatrixXd::Zero(dimW, p1a + p2a);

		Eigen::VectorXd alpha2 = Eigen::VectorXd::Zero(p1a + p2a);
		alpha2 << alpha21, alpha22; // diff bio, same alpha



		for (int i = 0; i < numSubj; i++) {
			//scalef = 0;
			Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);

			Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(p1a + p2a);
			Eigen::VectorXd bVec1 = Rcpp::as<Eigen::VectorXd>(bListElement[0]);
			Eigen::VectorXd bVec2 = Rcpp::as<Eigen::VectorXd>(bListElement[1]);


			bVeci << bVec1, bVec2;

			// Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
			Eigen::MatrixXd BAssociation = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
			Eigen::VectorXd latent = bVeci;
			double mu2, tausq;
			tausq = alpha2.transpose() * BAssociation * alpha2;

			Eigen::VectorXd w = W.row(i);
			//Eigen::VectorXd l = latent;
			Eigen::VectorXd l = BAssociation * alpha2 + bVeci;

			Eigen::MatrixXd wl = Eigen::MatrixXd::Zero(dimW, p1a + p2a);

			wl = w * l.transpose(); //

			mu2 = MultVV(w, gamma2) + MultVV(alpha2, bVeci);
			wwT = MultVVoutprod(w);
			//llT = MultVVoutprod(l);
			llT = MultVVoutprod(BAssociation * alpha2 + bVeci) + BAssociation;

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


			if (cmprsk(i) == 2)
			{

				scalefH02 = H02(risk2_index, 2);
				//SXX *= scalefH01;
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
			if (cmprsk(i) == 2) {
				Sw_inter += W.row(i);
				Sl_inter += latent;
				//Sl_inter +=  BAssociation  * alpha1 + bVeci;
			}
		}


		//NR update

		Sfull_inter = Eigen::VectorXd::Zero(dimW + p1a + p2a);
		Sfull_new = Eigen::VectorXd::Zero(dimW + p1a + p2a);
		info = Eigen::MatrixXd::Zero(dimW + p1a + p2a, dimW + p1a + p2a);

		Sfull_inter << Sw_inter, Sl_inter;
		//check2 = Sfull_inter;
		Sfull_new << Sw_new, Sl_new;

		// start row, start column, how many rows, how many col
		info.block(0, 0, dimW, dimW) = Sww_new;
		info.block(0, dimW, dimW, p1a + p2a) = Swl_new;
		info.block(dimW, 0, p1a + p2a, dimW) = Swl_new.transpose();
		info.block(dimW, dimW, p1a + p2a, p1a + p2a) = Sll_new;

		// NR update
		Eigen::VectorXd phi2 = Eigen::VectorXd::Zero(dimW + p1a + p2a);
		phi2 << gamma2, alpha2;

		phi2 += info.inverse() * (Sfull_inter - Sfull_new);



		Sw_new = Eigen::VectorXd::Zero(dimW);
		Sw_inter = Eigen::VectorXd::Zero(dimW);
		Sww_new = Eigen::MatrixXd::Zero(dimW, dimW);
		Sl_new = Eigen::VectorXd::Zero(p1a + p2a);
		Sl_inter = Eigen::VectorXd::Zero(p1a + p2a);
		Sll_new = Eigen::MatrixXd::Zero(p1a + p2a, p1a + p2a);
		Swl_new = Eigen::MatrixXd::Zero(dimW, p1a + p2a);


	return Rcpp::List::create(Rcpp::Named("FUNB") = FUNB,
		Rcpp::Named("beta1") = beta1New, Rcpp::Named("beta2") = beta2New, Rcpp::Named("beta") = betaNew,
		Rcpp::Named("Sig") = SigE,
		Rcpp::Named("sigma1") = sigma1, Rcpp::Named("sigma2") = sigma2,
		Rcpp::Named("H01") = H01, Rcpp::Named("H02") = H02,
		Rcpp::Named("phi1") = phi1, Rcpp::Named("phi2") = phi2);

}


