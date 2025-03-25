#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]


Rcpp::List getNoQuad(Rcpp::List XList, Rcpp::List YList, Rcpp::List ZList, Eigen::MatrixXd& W,
	Rcpp::List mdata, Rcpp::List mdataSList,
	Rcpp::List bList, Eigen::VectorXd sigmaInit, Rcpp::List sigmaiList,
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
	  
	  Eigen::MatrixXd XVXT = Eigen::MatrixXd::Zero(p, p);
	  Eigen::MatrixXd YZBX = Eigen::MatrixXd::Zero(p, 1);
	  Eigen::VectorXd betaNew = Eigen::VectorXd::Zero(p);
	  
	  for (int i = 0; i < numSubj; i++) {
	    
	    Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
	    Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
	    Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
	    Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
	    
	    Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[g]);
	    Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[g]);
	    Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[g]);
	    
	    double sigmag = sigmaInit(g);
	    
	    Eigen::MatrixXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
	    
	    XVXT = XVXT + Xtemp.transpose() * Xtemp / sigmag; //4x4
	    YZBX = YZBX + Xtemp.transpose() * (Ytemp - Ztemp * bVec) / sigmag; // 4x1
	    
	  }
	  
	  betaNew = XVXT.inverse() * YZBX;
	  betaFull.segment(index, p) = betaNew;
	  betaNewList[std::string("beta") + std::to_string(g+1)] = betaNew;
	  index += p;
	  
	  //~~~~~~~~~~~~
	  //
	  // sigma
	  //
	  // ~~~~~~~~~~~~~~~
	  
	  double numsig = 0;
	  int nijSum = 0;
	  

	  Eigen::MatrixXd ZZT = Eigen::MatrixXd::Zero(pRE, pRE);
	  pREindex = 0;
	  
	  for (int i = 0; i < numSubj; i++) {
	    Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
	    Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
	    Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
	    Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
	    
	    Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[g]);
	    Eigen::MatrixXd Ytemp = Rcpp::as<Eigen::MatrixXd>(yListElement[g]);
	    Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[g]);
	    
	    pRE = pREVec(g);
	    
	    Eigen::MatrixXd bVeci = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
	    Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
	    Eigen::MatrixXd sigig = sigmai.block(pREindex, pREindex, pRE, pRE);  // sigma i for gth biomarker
	    
	    Rcpp::List mdataList = Rcpp::as<Rcpp::List>(mdata[g]);
	    int numRep = Rcpp::as<int>(mdataList[i]);
	    
	    for (int nij = 0; nij < numRep; nij++) {

	      double epsilon = Ytemp(nij, 0) - MultVV(Xtemp.row(nij), betaNew);
	      double zb = MultVV(Ztemp.row(nij), bVeci);
	      
	      Eigen::MatrixXd ZZT = MultVVoutprod(Ztemp.row(nij));
	      Eigen::MatrixXd bbT = MultVVoutprod(bVeci);
	      
	      numsig += pow(epsilon, 2) - 2 * epsilon * zb + (ZZT * (sigig + bbT)).trace();
	   
	    }
	    
	    nijSum += numRep;
	  }
	  
	  sigmaVec(g) = numsig / nijSum;
	  pREindex += pRE;
	  
	  
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
	  
	  Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
	  Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
	  
	  for(int g = 0; g < numBio; g++){
	    Eigen::VectorXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
	    pRE = pREVec(g);
	    bVeci.segment(index, pRE) = bVec;
	    index += pRE;
	  }
	  
	  numSig += sigmai + MultVVoutprod(bVeci);
	  
	}
	
	SigE = numSig / numSubj;



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
		
		
		Eigen::VectorXd alpha1 = Eigen::VectorXd::Zero(pREtotal);
		Eigen::VectorXd alpha2 = Eigen::VectorXd::Zero(pREtotal);
		Rcpp::List alphaListElement1 = Rcpp::as<Rcpp::List>(alphaList[0]); // get first  risk
		Rcpp::List alphaListElement2 = Rcpp::as<Rcpp::List>(alphaList[1]); // get second risk
		
		index = 0;
		
		for(int g = 0; g < numBio; g++){
		  Eigen::VectorXd alpha1g = Rcpp::as<Eigen::VectorXd>(alphaListElement1[g]); // get alpha1// risk 1 bio g
		  Eigen::VectorXd alpha2g = Rcpp::as<Eigen::VectorXd>(alphaListElement2[g]); // get alpha2// risk 2 bio g
		  alpha1.segment(index, pRE) = alpha1g;
		  alpha2.segment(index, pRE) = alpha2g;
		  index += pRE;
		}
		
		index = 0;

		// // --- normal ----------------------
		for (int i = 0; i < numSubj; i++) {
			Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
				Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);


			for(int g = 0; g < numBio; g++){
			  Eigen::VectorXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
			  pRE = pREVec(g);
			  bVeci.segment(index, pRE) = bVec;
			  index += pRE;
			}
			
			index = 0;
			
			double muH1, tau1, tausq1;
			muH1 = MultVV(W.row(i), gamma1) + MultVV(alpha1, bVeci);
			tausq1 = alpha1.transpose() * sigmai * alpha1;

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
					  for(int g = 0; g < numBio; g++){
					    Eigen::VectorXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
					    pRE = pREVec(g);
					    bVeci.segment(index, pRE) = bVec;
					    index += pRE;
					  }
					  
					  index = 0;
					  
						sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
						muH1 = MultVV(W.row(i), gamma1) + MultVV(alpha1, bVeci);
						tausq1 = alpha1.transpose() * sigmai * alpha1;
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
		  index = 0;
		  
		  Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
		  Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
		  
		  for(int g = 0; g < numBio; g++){
		    Eigen::VectorXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
		    pRE = pREVec(g);
		    bVeci.segment(index, pRE) = bVec;
		    index += pRE;
		  }


		  index = 0;
		  
			double muH2, tausq2, tau2;
			muH2 = MultVV(W.row(i), gamma2) + MultVV(alpha2, bVeci);
			tausq2 = alpha2.transpose() * sigmai * alpha2;
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
					  for(int g = 0; g < numBio; g++){
					    Eigen::VectorXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
					    pRE = pREVec(g);
					    bVeci.segment(index, pRE) = bVec;
					    index += pRE;
					  }
					  
					  index = 0;
					  
					  
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
		Eigen::VectorXd Sl_new = Eigen::VectorXd::Zero(pREtotal);
		Eigen::VectorXd Sl_inter = Eigen::VectorXd::Zero(pREtotal);
		Eigen::MatrixXd Sll_new = Eigen::MatrixXd::Zero(pREtotal, pREtotal);
		Eigen::MatrixXd Swl_new = Eigen::MatrixXd::Zero(dimW, pREtotal);

		Eigen::VectorXd latent;

		Eigen::MatrixXd  wwT = Eigen::MatrixXd::Zero(dimW, dimW);
		Eigen::MatrixXd  llT = Eigen::MatrixXd::Zero(pREtotal, pREtotal);
		Eigen::MatrixXd  SwwT = Eigen::MatrixXd::Zero(dimW, dimW);
		Eigen::MatrixXd  SllT = Eigen::MatrixXd::Zero(pREtotal, pREtotal);
		Eigen::VectorXd  Sw = Eigen::MatrixXd::Zero(dimW, dimW);
		Eigen::VectorXd  Sl = Eigen::MatrixXd::Zero(pREtotal, pREtotal);
		Eigen::MatrixXd Swl = Eigen::MatrixXd::Zero(dimW, pREtotal);
		
		index = 0;

		for (int i = 0; i < numSubj; i++) {
			//scalef = 0;
			Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);

		  for(int g = 0; g < numBio; g++){
		    
		    Eigen::VectorXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
		    pRE = pREVec(g);
		    bVeci.segment(index, pRE) = bVec;
		    index += pRE;
		  }
		  
		  index = 0;

			// Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
			Eigen::MatrixXd BAssociation = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
			Eigen::VectorXd latent = bVeci;
			double mu1, tausq;
			tausq = alpha1.transpose() * BAssociation * alpha1;

			Eigen::VectorXd w = W.row(i);
			//Eigen::VectorXd l = latent;
			Eigen::VectorXd l = BAssociation * alpha1 + bVeci;

			Eigen::MatrixXd wl = Eigen::MatrixXd::Zero(dimW, pREtotal);

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
		
		index = 0;

		for (int i = 0; i < numSubj; i++)
		{
			// Eigen::VectorXd alpha1 = Eigen::VectorXd::Zero(pREtotal);
			// Eigen::VectorXd bVeci = Eigen::VectorXd::Zero(pREtotal);
			Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
			
			for(int g = 0; g < numBio; g++){
			  
			  Eigen::VectorXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
			  pRE = pREVec(g);
			  bVeci.segment(index, pRE) = bVec;
			  index += pRE;
			}
			
			index = 0;
			

			Eigen::VectorXd latent = bVeci;
			if (cmprsk(i) == 1) {
				Sw_inter += W.row(i);
				Sl_inter += latent;
			}
		}


		//NR update

		Eigen::VectorXd Sfull_inter = Eigen::VectorXd::Zero(dimW + pREtotal);
		Eigen::VectorXd Sfull_new = Eigen::VectorXd::Zero(dimW + pREtotal);
		Eigen::MatrixXd info = Eigen::MatrixXd::Zero(dimW + pREtotal, dimW + pREtotal);

		Sfull_inter << Sw_inter, Sl_inter;
		Sfull_new << Sw_new, Sl_new;

		// start row, start column, how many rows, how many col
		info.block(0, 0, dimW, dimW) = Sww_new;
		info.block(0, dimW, dimW, pREtotal) = Swl_new;
		info.block(dimW, 0, pREtotal, dimW) = Swl_new.transpose();
		info.block(dimW, dimW, pREtotal, pREtotal) = Sll_new;

		// NR update
		Eigen::VectorXd phi1 = Eigen::VectorXd::Zero(dimW + pREtotal);
		phi1 << gamma1, alpha1;
		phi1 += info.inverse() * (Sfull_inter - Sfull_new);


		Sw_new = Eigen::VectorXd::Zero(dimW);
		Sw_inter = Eigen::VectorXd::Zero(dimW);
		Sww_new = Eigen::MatrixXd::Zero(dimW, dimW);
		Sl_new = Eigen::VectorXd::Zero(pREtotal);
		Sl_inter = Eigen::VectorXd::Zero(pREtotal);
		Sll_new = Eigen::MatrixXd::Zero(pREtotal, pREtotal);
		Swl_new = Eigen::MatrixXd::Zero(dimW, pREtotal);
		// new matrix
		//Eigen::MatrixXd joint = Eigen::MatrixXd::Zero(p2 + p2, p2);
		// Eigen::MatrixXd  jjT = Eigen::MatrixXd::Zero(dimW   1a + p2a, dimW + pREtotal); // double check matrix dim might've switched
		// Eigen::MatrixXd SjjT = Eigen::MatrixXd::Zero(dimW + pREtotal, dimW + pREtotal);
		// Eigen::VectorXd Sj = Eigen::VectorXd::Zero(dimW + pREtotal);



		//     //MultVV(Z.row(j), FUNB.col(j) + tau * abscMat(0, j));

	//         //xb1 =  // fixed effect; contribute to LLH
	//         //zb = MultVV(FUNC, alpha1); // random effect component

		wwT = Eigen::MatrixXd::Zero(dimW, dimW);
		llT = Eigen::MatrixXd::Zero(pREtotal, pREtotal);
		SwwT = Eigen::MatrixXd::Zero(dimW, dimW);
		SllT = Eigen::MatrixXd::Zero(pREtotal, pREtotal);
		Sw = Eigen::MatrixXd::Zero(dimW, dimW);
		Sl = Eigen::MatrixXd::Zero(pREtotal, pREtotal);
		Swl = Eigen::MatrixXd::Zero(dimW, pREtotal);

		index = 0;


		for (int i = 0; i < numSubj; i++) {
			//scalef = 0;
			Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);

		  for(int g = 0; g < numBio; g++){
		    
		    Eigen::VectorXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
		    pRE = pREVec(g);
		    bVeci.segment(index, pRE) = bVec;
		    index += pRE;
		  }
		  
		  index = 0;

			// Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
			Eigen::MatrixXd BAssociation = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
			Eigen::VectorXd latent = bVeci;
			double mu2, tausq;
			tausq = alpha2.transpose() * BAssociation * alpha2;

			Eigen::VectorXd w = W.row(i);
			//Eigen::VectorXd l = latent;
			Eigen::VectorXd l = BAssociation * alpha2 + bVeci;

			Eigen::MatrixXd wl = Eigen::MatrixXd::Zero(dimW, pREtotal);

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

		index = 0;
		for (int i = 0; i < numSubj; i++)
		{

			Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
			for(int g = 0; g < numBio; g++){
			  
			  Eigen::VectorXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
			  pRE = pREVec(g);
			  bVeci.segment(index, pRE) = bVec;
			  index += pRE;
			}
			
			index = 0;


			Eigen::VectorXd latent = bVeci;
			if (cmprsk(i) == 2) {
				Sw_inter += W.row(i);
				Sl_inter += latent;
				//Sl_inter +=  BAssociation  * alpha1 + bVeci;
			}
		}


		//NR update

		Sfull_inter = Eigen::VectorXd::Zero(dimW + pREtotal);
		Sfull_new = Eigen::VectorXd::Zero(dimW + pREtotal);
		info = Eigen::MatrixXd::Zero(dimW + pREtotal, dimW + pREtotal);

		Sfull_inter << Sw_inter, Sl_inter;
		//check2 = Sfull_inter;
		Sfull_new << Sw_new, Sl_new;

		// start row, start column, how many rows, how many col
		info.block(0, 0, dimW, dimW) = Sww_new;
		info.block(0, dimW, dimW, pREtotal) = Swl_new;
		info.block(dimW, 0, pREtotal, dimW) = Swl_new.transpose();
		info.block(dimW, dimW, pREtotal, pREtotal) = Sll_new;

		// NR update
		Eigen::VectorXd phi2 = Eigen::VectorXd::Zero(dimW + pREtotal);
		phi2 << gamma2, alpha2;

		phi2 += info.inverse() * (Sfull_inter - Sfull_new);



		Sw_new = Eigen::VectorXd::Zero(dimW);
		Sw_inter = Eigen::VectorXd::Zero(dimW);
		Sww_new = Eigen::MatrixXd::Zero(dimW, dimW);
		Sl_new = Eigen::VectorXd::Zero(pREtotal);
		Sl_inter = Eigen::VectorXd::Zero(pREtotal);
		Sll_new = Eigen::MatrixXd::Zero(pREtotal, pREtotal);
		Swl_new = Eigen::MatrixXd::Zero(dimW, pREtotal);


		// Rcpp::Named("FUNB") = FUNB,
		
	return Rcpp::List::create(Rcpp::Named("beta") = betaFull,
                           Rcpp::Named("betaList") = betaNewList,
                           Rcpp::Named("sigmaVec") = sigmaVec,
		Rcpp::Named("Sig") = SigE,
		Rcpp::Named("H01") = H01, Rcpp::Named("H02") = H02,
		Rcpp::Named("phi1") = phi1, Rcpp::Named("phi2") = phi2);

}


