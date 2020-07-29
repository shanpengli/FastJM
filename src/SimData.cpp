//
//  SimData.cpp
//  OldJMmodel
//
//  Created by Shanpeng Li on 6/22/20.
//
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <string>
#include <iostream>
#include <fstream>

using namespace Rcpp;

void MulM(const gsl_matrix *XX, const gsl_vector *X, gsl_vector *beta)
{
    int p = XX->size1;
    int q = XX->size2;

    int i,j;
    double temp;

    for(i=0;i<p;i++)
    {
        temp=0;
        for(j=0;j<q;j++)  temp+=gsl_matrix_get(XX,i,j)*gsl_vector_get(X,j);
        gsl_vector_set(beta,i,temp);
    }

}

int GetN(double t)
{

    return (int)(t/1+1);

}

double Min(const double t1, const double t2)
{
    if(t1<t2) return t1;
    else return t2;
}

int gsl_linalg_cholesky_svx (const gsl_matrix * LLT,
                         gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      return 100;
    }
  else if (LLT->size2 != x->size)
    {
      return 100;
    }
  else
    {
      /* Solve for c using forward-substitution, L c = b */

      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);

      /* Perform back-substitution, U x = c */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LLT, x);

      return 0;
    }
}

int gsl_linalg_cholesky_solven (const gsl_matrix * LLT,
                           const gsl_vector * b,
                           gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      return 100;
    }
  else if (LLT->size1 != b->size)
    {
      return 100;
    }
  else if (LLT->size2 != x->size)
    {
      return 100;
    }
  else
    {
      /* Copy x <- b */

      gsl_vector_memcpy (x, b);

      /* Solve for c using forward-substitution, L c = b */

      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);

      /* Perform back-substitution, U x = c */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LLT, x);


      return 0;
    }
}



int gsl_linalg_cholesky_decompn (gsl_matrix * A)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      return 100;
    }
  else
    {
      size_t i,j,k;

      /* Do the first 2 rows explicitly.  It is simple, and faster.  And
       * one can return if the matrix has only 1 or 2 rows.
       */

      double A_00 = gsl_matrix_get (A, 0, 0);

      double L_00 = sqrt(A_00);

      if (A_00 <= 0)
        {
          return 100 ;
        }

      gsl_matrix_set (A, 0, 0, L_00);

      if (M > 1)
        {
          double A_10 = gsl_matrix_get (A, 1, 0);
          double A_11 = gsl_matrix_get (A, 1, 1);

          double L_10 = A_10 / L_00;
          double diag = A_11 - L_10 * L_10;
          double L_11 = sqrt(diag);

          if (diag <= 0)
            {
              return 100;
            }

          gsl_matrix_set (A, 1, 0, L_10);
          gsl_matrix_set (A, 1, 1, L_11);
        }

      for (k = 2; k < M; k++)
        {
          double A_kk = gsl_matrix_get (A, k, k);

          for (i = 0; i < k; i++)
            {
              double sum = 0;

              double A_ki = gsl_matrix_get (A, k, i);
              double A_ii = gsl_matrix_get (A, i, i);

              gsl_vector_view ci = gsl_matrix_row (A, i);
              gsl_vector_view ck = gsl_matrix_row (A, k);

              if (i > 0) {
                gsl_vector_view di = gsl_vector_subvector(&ci.vector, 0, i);
                gsl_vector_view dk = gsl_vector_subvector(&ck.vector, 0, i);

                gsl_blas_ddot (&di.vector, &dk.vector, &sum);
              }

              A_ki = (A_ki - sum) / A_ii;
              gsl_matrix_set (A, k, i, A_ki);
            }

          {
            gsl_vector_view ck = gsl_matrix_row (A, k);
            gsl_vector_view dk = gsl_vector_subvector (&ck.vector, 0, k);

            double sum = gsl_blas_dnrm2 (&dk.vector);
            double diag = A_kk - sum * sum;

            double L_kk = sqrt(diag);

            if (diag <= 0)
              {
                return 100;
              }

            gsl_matrix_set (A, k, k, L_kk);
          }
        }

      /* Now copy the transposed lower triangle to the upper triangle,
       * the diagonal is common.
       */

      for (i = 1; i < M; i++)
        {
          for (j = 0; j < i; j++)
            {
              double A_ij = gsl_matrix_get (A, i, j);
              gsl_matrix_set (A, j, i, A_ij);
            }
        }


    }

    return 0;
}

//' @export
// [[Rcpp::export]]

Rcpp::List SimData(SEXP k_val,SEXP p1_val,SEXP p1a_val,SEXP p2_val, SEXP g_val, SEXP truebeta, SEXP truegamma, SEXP truevee1, SEXP truevee2, SEXP randeffect, SEXP yfn, SEXP cfn, SEXP mfn)

{
    size_t k=as<int> (k_val);
    size_t p1=as<int> (p1_val);
    size_t p1a=as<int> (p1a_val);
    size_t p2=as<int> (p2_val);
    size_t g=as<int> (g_val);
    size_t n=0,n1,array_size=k*100;
    size_t i,j;

    std::string yfile=as<std::string>(yfn);
    std::string cfile=as<std::string>(cfn);
    std::string mfile=as<std::string>(mfn);
    const std::vector <double> tbeta_val=as<std::vector<double> >(truebeta);
    const std::vector <double> tgamma_val=as<std::vector<double> >(truegamma);
    const std::vector <double> reffect_val=as<std::vector<double> >(randeffect);
    const std::vector <double> tvee1_val=as<std::vector<double> >(truevee1);
    const std::vector <double> tvee2_val=as<std::vector<double> >(truevee2);
    if (tbeta_val.size() != p1)
    {
        Rprintf("Error in parameter p1 or truebeta.\n");
        return R_NilValue;
    }
    if (tgamma_val.size() != p2*2)
    {
        Rprintf("Error in parameter p2 (for two competing risks) or truegamma.\n");
        return R_NilValue;
    }

    gsl_matrix *C = gsl_matrix_calloc(k,p2+2);
    gsl_matrix *FY= gsl_matrix_calloc(array_size, p1+1);          /* initially assign a size for Y, call it fake Y ****/
    gsl_vector *M1= gsl_vector_calloc(k);
    gsl_matrix_set_zero(FY);

    /* allocate space for true parameters */
    gsl_vector *tbeta = gsl_vector_calloc(p1);
    gsl_vector *tvee1 = gsl_vector_calloc(p1a);
    gsl_vector *tvee2 = gsl_vector_calloc(p1a);
    gsl_matrix *tgamma = gsl_matrix_calloc(g,p2);

    gsl_matrix *VC= gsl_matrix_calloc(p1a,p1a);
    gsl_vector *RA= gsl_vector_calloc(p1a);
    gsl_vector *RI= gsl_vector_calloc(p1a);
    gsl_vector *S=gsl_vector_calloc(p1a);
    gsl_matrix *V=gsl_matrix_calloc(p1a,p1a);
    gsl_vector *W=gsl_vector_calloc(p1a);

    /* assign the true value to parameters */
    double tsigma = reffect_val[0];
    double tsigmab0 = reffect_val[1];
    double tsigmab1 = reffect_val[2];
    double rho = 0, sigmax = 0.1, prob = 0.5;
    double tsigmab01 = sqrt(tsigmab0*tsigmab1)*rho;
    for (size_t i=0;i<tbeta_val.size();i++) gsl_vector_set(tbeta,i,tbeta_val[i]);
    for (size_t i=0;i<g;i++)
    {
        for(size_t j=0;j<p2;j++)
        {
            gsl_matrix_set(tgamma,i,j,tgamma_val[(i)*(p2)+j]);
        }
    }

    for (size_t i=0;i<p1a;i++) gsl_vector_set(tvee1, i, tvee1_val[i]);
    for (size_t i=0;i<p1a;i++) gsl_vector_set(tvee2, i, tvee2_val[i]);

    gsl_matrix_set(VC,0,0,tsigmab0);
    gsl_matrix_set(VC,1,1,tsigmab1);
    gsl_matrix_set(VC,1,0,tsigmab01);

    for(i=0;i<p1a;i++)
    {
        for(j=i+1;j<p1a;j++)    gsl_matrix_set(VC,i,j,gsl_matrix_get(VC,j,i));
    }
    gsl_linalg_SV_decomp(VC,V,S,W);
    gsl_matrix_set_zero(V);
    for(i=0;i<p1a;i++) gsl_matrix_set(V,i,i,sqrt(gsl_vector_get(S,i)));

    int point=0;

    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    double lamda01=0.05, lamda02=0.1;               /* baseline hazard; constant over time */
    double temp,x1,x2,t1,t2,censor,u,b0,b1,error;
    double crate,rate1,rate2,max_censor=5;

    crate=0;
    rate1=0;
    rate2=0;

    for(j=0;j<k;j++)
        {
            for(i=0;i<p1a;i++)    gsl_vector_set(RI,i,gsl_ran_gaussian(r,1));
            MulM(V,RI,RA);
            MulM(VC,RA,RI);

            b0=gsl_vector_get(RI,0);
            b1=gsl_vector_get(RI,1);

            x1=gsl_ran_gaussian(r,sqrt(sigmax))+2;
            x2=gsl_ran_bernoulli(r,prob);
            censor=gsl_ran_exponential(r, 20);

            gsl_matrix_set(C,j,2,x1);
            gsl_matrix_set(C,j,3,x2);

            temp=gsl_matrix_get(tgamma,0,0)*x1+gsl_matrix_get(tgamma,0,1)*x2;
            temp=exp(temp+gsl_vector_get(tvee1,0)*b0+gsl_vector_get(tvee1,1)*b1);
            temp=lamda01*temp;
            t1=gsl_ran_exponential(r, 1/temp);

            temp=gsl_matrix_get(tgamma,1,0)*x1+gsl_matrix_get(tgamma,1,1)*x2;
            temp=exp(temp+gsl_vector_get(tvee2,0)*b0+gsl_vector_get(tvee2,1)*b1);
            temp=lamda02*temp;
            t2=gsl_ran_exponential(r, 1/temp);


            if(t1<=t2 && t1<=censor && t1<=max_censor)
            {
                rate1+=1;
                gsl_matrix_set(C,j,0,t1);
                gsl_matrix_set(C,j,1,1);
                n=GetN(t1);
            }

            if(t2<=t1 && t2<=censor && t2<=max_censor)
            {
                rate2+=1;
                gsl_matrix_set(C,j,0,t2);
                gsl_matrix_set(C,j,1,2);
                n=GetN(t2);
            }

            if(Min(censor, max_censor)<=t1 && Min(censor, max_censor)<=t2)
            {
                crate+=1;
                gsl_matrix_set(C,j,0,Min(max_censor,censor));
                gsl_matrix_set(C,j,1,0);
                n=GetN(Min(max_censor,censor));
            }


            for(i=point;i<point+n;i++)
            {
                gsl_matrix_set(FY,i,1,1.0);
                gsl_matrix_set(FY,i,3,x2);

                gsl_matrix_set(FY,i,2,(double)(i-point)*1);
                error=gsl_ran_gaussian(r,sqrt(tsigma));

                temp=gsl_vector_get(tbeta,0)+gsl_vector_get(tbeta,1)*gsl_matrix_get(FY,i,2)+gsl_vector_get(tbeta,2)*x2;
                temp=temp+error+b0+b1*gsl_matrix_get(FY,i,2);
                gsl_matrix_set(FY,i,0,temp);
            }

            point+=n;
            gsl_vector_set(M1,j,(double)n);

        }

        n1 = point;
        gsl_matrix *Y = gsl_matrix_calloc(n1,1+1+p1a+p1); /* true matrix for longitudinal outcome Y **/
        for(i=0;i<n1;i++)
        {
            gsl_matrix_set(Y,i,1,gsl_matrix_get(FY, i, 0));
            for (j=0;j<p1a;j++) gsl_matrix_set(Y,i,j+1+1,gsl_matrix_get(FY, i, j+1));
            for (j=0;j<p1;j++) gsl_matrix_set(Y, i, p1a+1+j+1, gsl_matrix_get(FY, i, j+1));
        }
        int q, p=0;

        for (j=0;j<k;j++)
        {
            u=(int)gsl_vector_get(M1,j);
            for (q=p;q<u+p;q++)
            {
                gsl_matrix_set(Y,q,0,j+1);
            }
            p=p+u;
        }
        /***** ###### output data ############## *******/

         //write Y file

         FILE *output_F;

         output_F=fopen(yfile.c_str(),"w");
         if (output_F == NULL)
         {
             Rprintf("Can't write y file\n");
         }
         fprintf(output_F,"%s   ", "ID response intercept_RE time_RE intercept time x1\n");


        for(i=0;i<Y->size1;i++)
        {
            for(j=0;j<Y->size2;j++)   fprintf(output_F,"%f   ", gsl_matrix_get(Y,i,j));
            if(i!=Y->size1-1)  fprintf(output_F,"\n");
        }
        //close writing files
        fclose(output_F);

        //write C file
        output_F=fopen(cfile.c_str(),"w");
        if (output_F == NULL)
        {
            Rprintf("Can't write C file\n");
            //return R_NilValue;
        }
        fprintf(output_F,"%s   ", "surv    failure_type   x1         x2\n");

        for(i=0;i<C->size1;i++)
        {
            for(j=0;j<C->size2;j++)   fprintf(output_F,"%f   ", gsl_matrix_get(C,i,j));
            if(i!=C->size1-1)  fprintf(output_F,"\n");
        }
        //close writing files
        fclose(output_F);

        //write M file

        output_F=fopen(mfile.c_str(),"w");
        if (output_F == NULL)
        {
            Rprintf("Can't write M file\n");
           // return R_NilValue;
        }

        for(i=0;i<M1->size;i++)
        {
            fprintf(output_F,"%f", gsl_vector_get(M1,i));
            if(i!=M1->size-1)  fprintf(output_F,"\n");
        }

        fclose(output_F);

        gsl_matrix_free(C);
        gsl_matrix_free(FY);
        gsl_matrix_free(Y);
        gsl_vector_free(M1);
        gsl_rng_free(r);

        gsl_matrix_free(VC);
        gsl_vector_free(RA);
        gsl_vector_free(RI);
        gsl_vector_free(S);
        gsl_matrix_free(V);
        gsl_vector_free(W);

        gsl_vector_free(tbeta);
        gsl_vector_free(tvee1);
        gsl_vector_free(tvee2);
        gsl_matrix_free(tgamma);

        Rcpp::List ret;
        ret["censoring_rate"] = crate/(double)k;
        ret["rate1"] = rate1/(double)k;
        ret["rate2"] = rate2/(double)k;
        ret["yfn"] = yfile.c_str();
        ret["cfn"] = cfile.c_str();
        ret["mfn"] = mfile.c_str();
        return ret;

}
