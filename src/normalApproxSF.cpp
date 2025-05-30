#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]


Rcpp::List normalApproxSF(Rcpp::List XList, Rcpp::List YList, Rcpp::List ZList, Eigen::MatrixXd& W,
	Rcpp::List mdata, Rcpp::List mdataSList,
	Rcpp::List bList, Eigen::VectorXd sigmaInit, Rcpp::List sigmaiList,
	Eigen::MatrixXd H01, Eigen::VectorXd& survtime, Eigen::VectorXd cmprsk,
	Eigen::VectorXd& gamma1, Rcpp::List alphaList,
	const Eigen::VectorXd& CUH01,
	const Eigen::VectorXd& HAZ01, const Eigen::MatrixXd& Sig,
	Rcpp::List betaList){


	Eigen::MatrixXd H01q = H01;

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

	    
	    double sigmag = sigmaInit(g);
	    
	    // Eigen::MatrixXd bVec = Rcpp::as<Eigen::VectorXd>(bListElement[g]);
	    
	    // std::cout << "bVecfull" << bVeci << std::endl;
	    // std::cout << "bVec" << bVecig << std::endl;
	    // std::cout << "pREindex" << pREindex << std::endl;
	    // std::cout << "pRE" << pRE << std::endl;
	    
	    XVXT = XVXT + Xtemp.transpose() * Xtemp / sigmag; //4x4
	    YZBX = YZBX + Xtemp.transpose() * (Ytemp - Ztemp * bVecig) / sigmag; // 4x1
	    
	  }
	  
	  betaNew = XVXT.inverse() * YZBX;
	  betaFull.segment(index, p) = betaNew;
	  // std::cout << "betaFull" <<betaFull << std::endl;
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
	  // pREindex = 0;
	  for (int i = 0; i < numSubj; i++) {
	    Rcpp::List xListElement = Rcpp::as<Rcpp::List>(XList[i]);
	    Rcpp::List yListElement = Rcpp::as<Rcpp::List>(YList[i]);
	    Rcpp::List zListElement = Rcpp::as<Rcpp::List>(ZList[i]);
	    // Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
	    
	    Eigen::MatrixXd Xtemp = Rcpp::as<Eigen::MatrixXd>(xListElement[g]);
	    Eigen::VectorXd Ytemp = Rcpp::as<Eigen::VectorXd>(yListElement[g]);
	    Eigen::MatrixXd Ztemp = Rcpp::as<Eigen::MatrixXd>(zListElement[g]);
	    
	    pRE = pREVec(g);
	    
	    Eigen::VectorXd bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
	    Eigen::VectorXd bVecig = bVeci.segment(pREindex, pRE);
	    
	    Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
	    Eigen::MatrixXd sigig = sigmai.block(pREindex, pREindex, pRE, pRE);  // sigma i for gth biomarker
	    
	    // std::cout << "bVeci" << bVeci << std::endl;
	    // std::cout << "bVecig" << bVecig << std::endl;
	    // std::cout << "sigmai" << sigmai << std::endl;
	    // std::cout << "sigig" << sigig<< std::endl;
	    
	    Rcpp::List mdataList = Rcpp::as<Rcpp::List>(mdata[g]);
	    int numRep = Rcpp::as<int>(mdataList[i]);
	    
	    for (int nij = 0; nij < numRep; nij++) {

	      double epsilon = Ytemp(nij) - MultVV(Xtemp.row(nij), betaNew);
	      double zb = MultVV(Ztemp.row(nij), bVecig);
	      
	      Eigen::MatrixXd ZZT = MultVVoutprod(Ztemp.row(nij));
	      Eigen::MatrixXd bbT = MultVVoutprod(bVecig);
	      
	      
	      // std::cout << "zzt" << ZZT << std::endl;
	      // std::cout << "bbt" << bbT << std:: endl;

	      numsig += pow(epsilon, 2) - 2 * epsilon * zb + (ZZT * (sigig + bbT)).trace();
	      
	      // if(g == 1){
	        // std::cout << "z" << Ztemp.row(nij) << std::endl;
	        // std::cout << "bVecig" << bVecig << std::endl;
	        // std::cout << "zb" << zb << std::endl;
	        // std::cout << "zzt" << ZZT << std::endl;
	        // std::cout << "bbt" << bbT << std:: endl;
	   //    std::cout << "epsilon" << epsilon << std::endl;
	   //    std::cout << "middle part" << 2 * epsilon * zb << std:: endl;
	   // std::cout << "trace thing" << (ZZT * (sigig + bbT)).trace() << std::endl;
	   // // std::cout << "z" << Ztemp.row(nij) << std::endl;
	   // std::cout << "bVecig" << bVecig << std::endl;
	   // std::cout << "numsig" <<numsig << std::endl;
	      // }
	   
	    }
	    
	    nijSum += numRep;
	  }
	  
	  sigmaVec(g) = numsig / nijSum;
	  pREindex += pRE;
	  
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
	  
	  // Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);
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
		double denomH = 0;

		double dem1 = 0;
		double dem2 = 0;
		int a = H01.rows();
		int risk1_index = a - 1;
		
		Eigen::VectorXd alpha1 = Eigen::VectorXd::Zero(pREtotal);
		Rcpp::List alphaListElement1 = Rcpp::as<Rcpp::List>(alphaList[0]); // get first  risk
		
		index = 0;
		
		for(int g = 0; g < numBio; g++){
		  pRE = pREVec(g);
		  Eigen::VectorXd alpha1g = Rcpp::as<Eigen::VectorXd>(alphaListElement1[g]); // get alpha1// risk 1 bio g
		  alpha1.segment(index, pRE) = alpha1g;
		  // std::cout << "alpha 1g " << alpha1g<< std::endl;
		  // std::cout << "alpha 2g " << alpha2g<< std::endl;
		  // std::cout << "alpha 1 " << alpha1<< std::endl;
		  // std::cout << "alpha 2 " << alpha2<< std::endl;
		  // std::cout<< "index " << index << std::endl;
		  // std::cout<< "pRE " << pRE << std::endl;
		  index += pRE;

		}
		
		index = 0;

		// // --- normal ----------------------
		
		
		// HAZARD 1
		for (int i = 0; i < numSubj; i++) {
				Eigen::MatrixXd sigmai = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
			
			  Eigen::VectorXd bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
			// Eigen::VectorXd bVecig = bVeci.segment(pREindex, pRE);
			
			index = 0;
			
			double muH1, tausq1;
			muH1 = MultVV(W.row(i), gamma1) + MultVV(alpha1, bVeci);
			tausq1 = alpha1.transpose() * sigmai * alpha1;

			dem1 += exp(muH1 + 0.5 * tausq1);
			
			// std::cout << "muH1" << mu << std::endl;
			// std::cout << "dem1" << dem1 << std::endl;
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
					
					  bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
					  
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

		/////////////////////////////////////////
		
		/////////////////////////////////////////
		
		// PHI
		
		////

		double scalefH01 = 0;
		double scalef;

		risk1_index = a - 1;

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
			// Rcpp::List bListElement = Rcpp::as<Rcpp::List>(bList[i]);

		  // for(int g = 0; g < numBio; g++){
		  //   Eigen::VectorXd bVec = Rcpp::as<Eigen::VectorXd>(bList[i]);
		  //   pRE = pREVec(g);
		  //   index += pRE;
		  // }
		  
		  Eigen::VectorXd bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
		  
		  // index = 0;

			Eigen::MatrixXd BAssociation = Rcpp::as<Eigen::MatrixXd>(sigmaiList[i]);
			Eigen::VectorXd latent = bVeci;
			double mu1, tausq;
			tausq = alpha1.transpose() * BAssociation * alpha1;

			Eigen::VectorXd w = W.row(i);
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
				// 
				//         SwwT /= scalefH01;
				//         SllT /= scalefH01;
				//         Swl /= scalefH01;
				// 
				//          // s
				//        Sw *= scalefH01;
				//        Sl *= scalefH01;
				//         Sw_new += Sw;
				//         Sl_new += Sl;
				//         Sw /= scalefH01;
				//         Sl /= scalefH01;
				// 
				// 
				//         risk1_index--;
				// 
				//     }
				// 
				//     else if (survtime(i + 1) != survtime(i))
				//     {
				//         scalefH01 = H01(risk1_index, 2);
				// 
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
				// 
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
				//         // for (i = i + 1; i < numSubj; i++)
				//         // {
				//         // 
				//         //     bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
				//         //     latent = bVeci;
				//         //     w = W.row(i);
				//         //     l = latent;
				//         // 
				//         // 
				//         //     mu1 = MultVV(w, gamma1) + MultVV(alpha1, l);
				//         //     wwT =  MultVVoutprod(W.row(i));
				//         //     llT = MultVVoutprod(latent);
				// 
				// //             scalef = exp(mu1);
				// 
				// //             wl = Eigen::MatrixXd::Zero(dimW, p1a+p2a);
				// 
				// //             wl = w * l.transpose();
				// 
				// // // for I
				// // wwT *= scalef;
				// // llT *= scalef;
				// // wl *= scalef;
				// // SwwT += wwT;
				// // SllT += llT;
				// // Swl += wl;
				// 
				// 
				// // // // for S
				// // // wwT *= scalef;
				// // // llT *= scalef;
				// // // SwwT += wwT;
				// // // SllT += llT;
				// 
				// 
				// // // for S
				// // w *= scalef;
				// // l *= scalef;
				// // Sw += w;
				// // Sl += l;
				// 
				// 
				// //             if (i == numSubj - 1)
				// //             {
				// //                 scalefH01 = H01(risk1_index, 2);
				// //         SwwT *= scalefH01;
				// //         SllT *= scalefH01;
				// //         Swl *= scalefH01;
				// //         Sww_new += SwwT;
				// //         Sll_new += SllT;
				// //         Swl_new += Swl;
				// //         SwwT /= scalefH01;
				// //         SllT /= scalefH01;
				// //         Swl /= scalefH01;
				// 
				// //          // s
				// //        Sw *= scalefH01;
				// //        Sl *= scalefH01;
				// //         Sw_new += Sw;
				// //         Sl_new += Sl;
				// //         Sw /= scalefH01;
				// //         Sl /= scalefH01;
				// 
				// //         risk1_index--;
				// //                 break;
				// //             }
				// //             else if (survtime(i + 1) != survtime(i))
				// //             {
				// //                 scalefH01 = H01(risk1_index, 2);
				// 
				// //                 // i
				// //         SwwT *= scalefH01;
				// //         SllT *= scalefH01;
				// //         Swl *= scalefH01;
				// //         Sww_new += SwwT;
				// //         Sll_new += SllT;
				// //         Swl_new += Swl;
				// //         SwwT /= scalefH01;
				// //         SllT /= scalefH01;
				// //         Swl /= scalefH01;
				// 
				// //          // s
				// //        Sw *= scalefH01;
				// //        Sl *= scalefH01;
				// //         Sw_new += Sw;
				// //         Sl_new += Sl;
				// //         Sw /= scalefH01;
				// //         Sl /= scalefH01;
				// //         risk1_index--;
				// //                 break;
				// //             }
				// //             else continue;
				// //         }
				// //     }
				// }
				// // else continue;
			}

		}
		
		index = 0;

		for (int i = 0; i < numSubj; i++)
		{
			
			Eigen::VectorXd bVeci = Rcpp::as<Eigen::VectorXd>(bList[i]);
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
		// std::cout<< "gamma1 " <<  gamma1 << std::endl;
		// std::cout<< "info " << info.inverse() << std::endl;
		phi1 << gamma1, alpha1;
		phi1 += info.inverse() * (Sfull_inter - Sfull_new);

		// Rcpp::Named("FUNB") = FUNB,
		
	return Rcpp::List::create(Rcpp::Named("beta") = betaFull,
                           Rcpp::Named("betaList") = betaNewList,
                           Rcpp::Named("sigmaVec") = sigmaVec,
		                        Rcpp::Named("Sig") = SigE,
		                    Rcpp::Named("H01") = H01,
		                    Rcpp::Named("phi1") = phi1);

}


