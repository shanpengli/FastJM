//
//  jmcsf.cpp
//  FastJM
//
//  Created by Shanpeng Li on 6/16/20.
//

#include "jmcsf.hpp"
#include "comp.h"
namespace jmcsfspace {

    double GetPosbi(
                    const gsl_matrix *Y,
                    const gsl_vector *beta,
                    const gsl_vector *M1,
                    const gsl_matrix *Sigcov,
                    const double sigma,
                    gsl_matrix *Xinv,
                    gsl_matrix *Posbi,
                    const int k,
                    const int p1,
                    const int p1a
                    )
    {
        int i=0,f=0,w=0,t,j;
        for (j=0;j<k;j++)
        {
            i=(int)gsl_vector_get(M1,j);
            gsl_matrix *Z = gsl_matrix_calloc(i, p1a);
            gsl_matrix *ZT = gsl_matrix_calloc(p1a, i);
            gsl_matrix *ZZT = gsl_matrix_calloc(i, i);
            gsl_matrix *DZT = gsl_matrix_calloc(p1a, i);
            gsl_matrix *DZT1 = gsl_matrix_calloc(p1a, i);
            gsl_matrix *V = gsl_matrix_calloc(i, i);
            gsl_matrix *X =gsl_matrix_calloc(i, p1);
            gsl_matrix *XT = gsl_matrix_calloc(p1, i);
            gsl_matrix *XTV = gsl_matrix_calloc(p1, i);
            gsl_matrix *XTX = gsl_matrix_calloc(p1, p1);
            gsl_matrix_set_identity(V);
            gsl_vector *res = gsl_vector_calloc(i);
            gsl_vector *y = gsl_vector_calloc(i);
            gsl_vector *bi = gsl_vector_calloc(p1a);
            for(w=f;w<i+f;w++)
            {
                for (t=0;t<p1a;t++) gsl_matrix_set(Z, w-f, t, gsl_matrix_get(Y, w, t+1));
                for (t=0;t<p1;t++) gsl_matrix_set(X, w-f, t, gsl_matrix_get(Y, w, t+1+p1a));
                gsl_vector_set(y, w-f, gsl_matrix_get(Y, w, 0));
            }

            TransM(Z, ZT);
            TransM(X, XT);
            MulM(X, beta, res);
            gsl_vector_sub(y, res);
            MulMM(Sigcov, ZT, DZT);
            MulMM(Z, DZT, ZZT);
            gsl_matrix_scale(V, sigma);
            gsl_matrix_add(V, ZZT);
            inv_matrix(V);
            MulMM(DZT, V, DZT1);
            MulM(DZT1, y, bi);
            for (t=0;t<p1a;t++) gsl_matrix_set(Posbi, j, t, gsl_vector_get(bi, t));
            MulMM(XT, V, XTV);
            MulMM(XTV, X, XTX);
            gsl_matrix_add(Xinv, XTX);
            f=f+i;
            gsl_matrix_free(Z);
            gsl_matrix_free(ZT);
            gsl_matrix_free(ZZT);
            gsl_matrix_free(DZT);
            gsl_matrix_free(DZT1);
            gsl_matrix_free(V);
            gsl_matrix_free(X);
            gsl_matrix_free(XT);
            gsl_matrix_free(XTV);
            gsl_matrix_free(XTX);
            gsl_vector_free(res);
            gsl_vector_free(y);
            gsl_vector_free(bi);
        }
        return 0;
    }

    double GetPoscov(
                     const gsl_matrix *Y,
                     const gsl_vector *beta,
                     const gsl_vector *M1,
                     const gsl_matrix *Sigcov,
                     const double sigma,
                     const gsl_matrix *Xinv,
                     gsl_matrix *Poscov,
                     const int k,
                     const int p1,
                     const int p1a
                     )
    {
        int i=0;
        int f=0,w=0,t,j,u;
        gsl_permutation *pp = gsl_permutation_calloc(p1a);
         for (j=0;j<k;j++)
         {
             i=(int)gsl_vector_get(M1,j);
             gsl_matrix *Z = gsl_matrix_calloc(i, p1a);
             gsl_matrix *ZT = gsl_matrix_calloc(p1a, i);
             gsl_matrix *ZZT = gsl_matrix_calloc(i, i);
             gsl_matrix *DZT = gsl_matrix_calloc(p1a, i);
             gsl_matrix *ZDT = gsl_matrix_calloc(i, p1a);
             gsl_matrix *V = gsl_matrix_calloc(i, i);
             gsl_matrix *X =gsl_matrix_calloc(i, p1);
             gsl_matrix *VX = gsl_matrix_calloc(i, p1);
             gsl_matrix *XT = gsl_matrix_calloc(p1, i);
             gsl_matrix *XXT = gsl_matrix_calloc(i, i);
             gsl_matrix *XXTV = gsl_matrix_calloc(i, i);
             gsl_matrix *cov = gsl_matrix_calloc(p1a, p1a);
             gsl_matrix *covT = gsl_matrix_calloc(p1a, p1a);
             gsl_matrix_set_identity(V);
             gsl_vector *bi = gsl_vector_calloc(p1a);
             for(w=f;w<i+f;w++)
             {
                 for (t=0;t<p1a;t++) gsl_matrix_set(Z, w-f, t, gsl_matrix_get(Y, w, t+1));
                 for (t=0;t<p1;t++) gsl_matrix_set(X, w-f, t, gsl_matrix_get(Y, w, t+1+p1a));
             }
             TransM(Z, ZT);
             TransM(X, XT);
             MulMM(Sigcov, ZT, DZT);
             MulMM(Z, DZT, ZZT);
             gsl_matrix_scale(V, sigma);
             gsl_matrix_add(V, ZZT);
             inv_matrix(V);
             MulMM(V, X, VX);
             MulMM(VX, Xinv, X);
             MulMM(X, XT, XXT);
             MulMM(XXT, V, XXTV);
             gsl_matrix_sub(V, XXTV);
             MulMM(DZT, V, ZT);
             TransM(DZT, ZDT);
             MulMM(ZT, ZDT, cov);
             gsl_matrix_scale(cov, -1);
             gsl_matrix_add(cov, Sigcov);
             inv_matrix(cov);
             gsl_linalg_cholesky_decomp(cov);
             if (p1a>1)
             {
                for(u=1;u<p1a;u++)
                 {
                     for(t=0;t<p1a-u;t++)   gsl_matrix_set(cov,t,u+t,0);
                 }
             }
             TransM(cov, covT);
             gsl_linalg_LU_invert(covT, pp, cov);
             for (t=0;t<p1a;t++)
             {
                 gsl_matrix_get_row(bi, cov, t);
                 gsl_matrix_set_row(Poscov, j*p1a+t, bi);
             }
             f=f+i;
             gsl_matrix_free(Z);
             gsl_matrix_free(ZT);
             gsl_matrix_free(ZZT);
             gsl_matrix_free(DZT);
             gsl_matrix_free(ZDT);
             gsl_matrix_free(V);
             gsl_matrix_free(X);
             gsl_matrix_free(VX);
             gsl_matrix_free(XT);
             gsl_matrix_free(XXT);
             gsl_matrix_free(XXTV);
             gsl_matrix_free(cov);
             gsl_matrix_free(covT);
             gsl_vector_free(bi);
         }
        gsl_permutation_free(pp);
        return 0;
    }

    int EM(
           gsl_vector *beta,
           gsl_matrix *gamma,
           gsl_vector *vee1,
           gsl_matrix *H01,
           double *sigma,
           gsl_matrix *sig,
           const int p1a,
           const gsl_matrix *Y,
           const gsl_matrix *C,
           const gsl_vector *M1,
           const gsl_matrix *Posbi,
           const gsl_matrix *Poscov,
           const int point,
           const std::vector<double> xs,
           const std::vector<double> ws
           )

    {

        int p1=beta->size;
        int p2=gamma->size2;
        int g =gamma->size1;
        int a =H01->size2;

        int n1 = Y->size1;
        int k = M1->size;

        int p,q,j,t,u;

        gsl_matrix *FUNB=gsl_matrix_calloc(p1a,k),
                   *FUNBS=gsl_matrix_calloc(p1a*(p1a+1)/2,k),
                   *FUNBSE=gsl_matrix_calloc(g*p1a*(p1a+1)/2,k),
                   *FUNBE=gsl_matrix_calloc(g*p1a,k),
                   *FUNE=gsl_matrix_calloc(g,k);

        int status;
        status = GetE(FUNB,FUNBS,FUNBSE,FUNBE,FUNE,beta,gamma,vee1,H01,*sigma,sig,Y,C,M1,Posbi,Poscov,p1a,point,xs,ws);
        if (status==100) return status;

        gsl_vector * SX = gsl_vector_calloc(p2);
        gsl_matrix * SXX = gsl_matrix_calloc(p2,p2);
        gsl_matrix * XX = gsl_matrix_calloc(p2,p2);
        gsl_vector * X = gsl_vector_calloc(p2);
        gsl_vector * gammai = gsl_vector_calloc(p2);

        double scalef;

        gsl_vector * Z = gsl_vector_calloc(p1);                         /* covariates for Y */
        gsl_vector * SZ = gsl_vector_calloc(p1);
        gsl_matrix * SZZ = gsl_matrix_calloc(p1,p1);                    /* store sum of ZZ' */
        gsl_matrix * ZZ = gsl_matrix_calloc(p1,p1);                      /* store ZZ' */

        gsl_vector * bi = gsl_vector_calloc(p1a);
        gsl_matrix * bs = gsl_matrix_calloc(p1a,p1a);
        gsl_vector * xtilde = gsl_vector_calloc(p1a);
        gsl_vector * xtilde1 = gsl_vector_calloc(p1a);

        gsl_matrix *D=gsl_matrix_calloc(p1a,p1a);
        gsl_vector *N=gsl_vector_calloc(p1a);

        gsl_matrix *TD=gsl_matrix_calloc(p1a,p1a);
        gsl_vector *TN=gsl_vector_calloc(p1a);
        gsl_matrix *TDD=gsl_matrix_calloc(p1a,p1a);
        gsl_vector *TNN=gsl_vector_calloc(p1a);

        gsl_vector * SX_inter = gsl_vector_calloc(p2);
        gsl_vector * SX_new = gsl_vector_calloc(p2);
        gsl_matrix * SXX_new = gsl_matrix_calloc(p2,p2);
        gsl_vector * X_new = gsl_vector_calloc(p2);

        /* calculate beta, sig */

        gsl_matrix_set_zero(SZZ);
        gsl_vector_set_zero(SZ);
        gsl_matrix_set_zero(sig);

        p=0;

        for(j=0;j<k;j++)
        {
            u=(int)gsl_vector_get(M1,j);

            for(q=p;q<u+p;q++)
            {
                for(t=0;t<p1;t++)   gsl_vector_set(Z,t,gsl_matrix_get(Y,q,1+p1a+t));
                MulV(Z,ZZ);
                gsl_matrix_add(SZZ,ZZ);

                for(t=0;t<p1a;t++)   gsl_vector_set(xtilde,t,gsl_matrix_get(Y,q,1+t));
                for(t=0;t<p1a;t++)   gsl_vector_set(bi,t,gsl_matrix_get(FUNB,t,j));
                gsl_vector_scale(Z,gsl_matrix_get(Y,q,0)-MulVV(bi,xtilde));
                gsl_vector_add (SZ,Z);
            }

            p=p+u;

            for(t=0;t<p1a;t++)  gsl_matrix_set(sig,t,t,gsl_matrix_get(sig,t,t)+gsl_matrix_get(FUNBS,t,j));
            if (p1a>1)
            {
                for(q=1;q<p1a;q++)
                {
                    for(t=0;t<p1a-q;t++)   gsl_matrix_set(sig,t,q+t,gsl_matrix_get(sig,t,q+t)+gsl_matrix_get(FUNBS,p1a+t+(q-1)*(p1a-1),j));
                }
            }

        }

        gsl_matrix_scale(sig,1/(double)k);

        status=inv_matrix(SZZ);
        if(status==100)    return status;
        MulM(SZZ,SZ,beta);

        /* calculate sigma */
        *sigma=0;
        p=0;
        int i;
        for(j=0;j<k;j++)
        {
            u=(int)gsl_vector_get(M1,j);

            for(q=p;q<u+p;q++)
            {
                for(t=0;t<p1;t++)   gsl_vector_set(Z,t,gsl_matrix_get(Y,q,1+p1a+t));
                for(t=0;t<p1a;t++)   gsl_vector_set(bi,t,gsl_matrix_get(FUNB,t,j));
                for(t=0;t<p1a;t++)   gsl_vector_set(xtilde,t,gsl_matrix_get(Y,q,1+t));
                for(t=0;t<p1a;t++)   gsl_matrix_set(bs,t,t,gsl_matrix_get(FUNBS,t,j));

                if (p1a>1)
                {
                    for(i=1;i<p1a;i++)
                    {
                        for(t=0;t<p1a-i;t++)   gsl_matrix_set(bs,t,i+t,gsl_matrix_get(FUNBS,p1a+t+(i-1)*(p1a-1),j));
                    }
                    for(t=0;t<p1a;t++)
                    {
                        for(i=0;i<t;i++)   gsl_matrix_set(bs,t,i,gsl_matrix_get(bs,i,t));
                    }
                }
                MulM(bs,xtilde,xtilde1);
            *sigma+=gsl_pow_2(gsl_matrix_get(Y,q,0)-MulVV(Z,beta))-2*(gsl_matrix_get(Y,q,0)-MulVV(Z,beta))
                *MulVV(bi,xtilde)+MulVV(xtilde,xtilde1);
            }
            p=p+u;
        }

        *sigma=*sigma/(double)n1;



        /* calculate H01*/


        double dem;
        int risk1_index=a-1;
        /*New risk set calculation*/
         dem=0;
            for (j=0;j<k;j++)
            {
                for(u=0;u<p2;u++)
                {
                    gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
                    gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));
                }
                dem = dem + gsl_matrix_get(FUNE,0,j)*exp(MulVV(X,gammai));

                if (gsl_matrix_get(C,j,1) == 1)
                {

                    if (j == k-1)
                    {
                        //dem_1[risk1_index] = dem;
                        gsl_matrix_set(H01,2,risk1_index,gsl_matrix_get(H01,1,risk1_index)/dem);
                        risk1_index--;
                    }
                    else if (gsl_matrix_get(C,j+1,0) != gsl_matrix_get(C,j,0))
                    {
                        //dem_1[risk1_index] = dem;
                        gsl_matrix_set(H01,2,risk1_index,gsl_matrix_get(H01,1,risk1_index)/dem);
                        risk1_index--;
                    }

                    else
                    {
                        for (j=j+1;j<k;j++)
                        {
                            for(u=0;u<p2;u++)
                            {
                                gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
                                gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));
                            }
                            dem = dem + gsl_matrix_get(FUNE,0,j)*exp(MulVV(X,gammai));
                            if (j == k-1)
                            {
                                gsl_matrix_set(H01,2,risk1_index,gsl_matrix_get(H01,1,risk1_index)/dem);
                                risk1_index--;
                                break;
                            }
                            else if (gsl_matrix_get(C,j+1,0) != gsl_matrix_get(C,j,0))
                            {
                                //dem_1[risk1_index] = dem;
                                gsl_matrix_set(H01,2,risk1_index,gsl_matrix_get(H01,1,risk1_index)/dem);
                                risk1_index--;
                                break;
                            }
                            else continue;
                        }
                    }

                }
             else continue;
            }


        /* calculate gamma */
        /*New algorithm*/
        double scalefH01;
        risk1_index=a-1;
        gsl_matrix_set_zero(SXX);
        gsl_vector_set_zero(SX);
        gsl_matrix_set_zero(SXX_new);
        gsl_vector_set_zero(SX_new);
        gsl_vector_set_zero(SX_inter);
        for (j=0; j<k; j++)
        {
           for(u=0;u<p2;u++)
           {
               gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
               gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));
           }
           MulV(X,XX);
           scalef = gsl_matrix_get(FUNE,0,j)*exp(MulVV(X,gammai));

           gsl_matrix_scale(XX, scalef);
           gsl_matrix_add(SXX, XX);
           gsl_vector_scale(X, scalef);
           gsl_vector_add(SX, X);

           if (gsl_matrix_get(C,j,1) == 1)
           {
               if (j == k-1)
               {
                   scalefH01 = gsl_matrix_get(H01, 2, risk1_index);
                   gsl_matrix_scale(SXX, scalefH01);
                   gsl_matrix_add(SXX_new, SXX);
                   gsl_matrix_scale(SXX, 1/scalefH01);
                   gsl_vector_scale(SX, scalefH01);
                   gsl_vector_add(SX_new, SX);
                   gsl_vector_scale(SX, 1/scalefH01);
                   risk1_index--;
               }
               else if (gsl_matrix_get(C,j+1,0) != gsl_matrix_get(C,j,0))
               {
                   scalefH01 = gsl_matrix_get(H01, 2, risk1_index);
                   gsl_matrix_scale(SXX, scalefH01);
                   gsl_matrix_add(SXX_new, SXX);
                   gsl_matrix_scale(SXX, 1/scalefH01);
                   gsl_vector_scale(SX, scalefH01);
                   gsl_vector_add(SX_new, SX);
                   gsl_vector_scale(SX, 1/scalefH01);
                   risk1_index--;
               }
               else
               {
                   for (j=j+1;j<k;j++)
                   {
                       for(u=0;u<p2;u++)
                       {
                           gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
                           gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));
                       }
                       MulV(X,XX);
                       scalef = gsl_matrix_get(FUNE,0,j)*exp(MulVV(X,gammai));
                       gsl_matrix_scale(XX, scalef);
                       gsl_matrix_add(SXX, XX);
                       gsl_vector_scale(X, scalef);
                       gsl_vector_add(SX, X);

                       if (j == k-1)
                       {
                           scalefH01 = gsl_matrix_get(H01, 2, risk1_index);
                           gsl_matrix_scale(SXX, scalefH01);
                           gsl_matrix_add(SXX_new, SXX);
                           gsl_matrix_scale(SXX, 1/scalefH01);
                           gsl_vector_scale(SX, scalefH01);
                           gsl_vector_add(SX_new, SX);
                           gsl_vector_scale(SX, 1/scalefH01);
                           risk1_index--;
                           break;
                       }
                       else if (gsl_matrix_get(C,j+1,0) != gsl_matrix_get(C,j,0))
                       {
                           scalefH01 = gsl_matrix_get(H01, 2, risk1_index);
                           gsl_matrix_scale(SXX, scalefH01);
                           gsl_matrix_add(SXX_new, SXX);
                           gsl_matrix_scale(SXX, 1/scalefH01);
                           gsl_vector_scale(SX, scalefH01);
                           gsl_vector_add(SX_new, SX);
                           gsl_vector_scale(SX, 1/scalefH01);
                           risk1_index--;
                           break;
                       }
                       else continue;
                   }
               }
           }
           else continue;
        }

        for (j=0;j<k;j++)
        {
            for(u=0;u<p2;u++)
            {
                gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
            }
            if ((int)gsl_matrix_get(C,j,1) == 1)
            {
                gsl_vector_add(SX_inter, X);
            }

        }
        gsl_vector_scale(SX_new, -1);
        gsl_vector_add(SX_inter, SX_new);
        status=inv_matrix(SXX_new);
        if(status==100)  return status;
        MulM(SXX_new,SX_inter,X_new);

        for(j=0;j<p2;j++) gsl_matrix_set(gamma,0,j,gsl_matrix_get(gamma,0,j)+gsl_vector_get(X_new,j));

        /* calculate vee */
        gsl_matrix_set_zero(TD);
        gsl_vector_set_zero(TN);
        risk1_index = a-1;

        for (j=0;j<k;j++)
        {
            for(t=0;t<p1a;t++)   gsl_matrix_set(D,t,t,gsl_matrix_get(FUNBSE,t,j));

            if(p1a>1)
            {
                for(i=1;i<p1a;i++)
                {
                    for(t=0;t<p1a-i;t++)   gsl_matrix_set(D,t,i+t,gsl_matrix_get(FUNBSE,p1a+t+(i-1)*(p1a-1),j));
                }
            }
            for(t=0;t<p1a;t++)
            {
                for(i=0;i<t;i++)   gsl_matrix_set(D,t,i,gsl_matrix_get(D,i,t));
            }

            for (t=0;t<p1a;t++) gsl_vector_set(N,t,gsl_matrix_get(FUNBE,t,j));

            for(u=0;u<p2;u++)
            {
                gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
                gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));
            }
            gsl_matrix_scale(D,exp(MulVV(X,gammai)));
            gsl_matrix_add(TD,D);
            gsl_vector_scale(N,exp(MulVV(X,gammai)));
            gsl_vector_add(TN,N);

            if (gsl_matrix_get(C,j,1) == 1)
            {
              if (j == k-1)
              {
                gsl_matrix_scale(TD, gsl_matrix_get(H01, 2, risk1_index));
                gsl_matrix_add(TDD,TD);
                gsl_matrix_scale(TD, 1/gsl_matrix_get(H01, 2, risk1_index));
                gsl_vector_scale(TN, gsl_matrix_get(H01, 2, risk1_index));
                gsl_vector_add(TNN,TN);
                gsl_vector_scale(TN, 1/gsl_matrix_get(H01, 2, risk1_index));
                risk1_index--;


              }
              else if (gsl_matrix_get(C,j+1,0) != gsl_matrix_get(C,j,0))
              {
                 gsl_matrix_scale(TD, gsl_matrix_get(H01, 2, risk1_index));
                 gsl_matrix_add(TDD,TD);
                 gsl_matrix_scale(TD, 1/gsl_matrix_get(H01, 2, risk1_index));
                 gsl_vector_scale(TN, gsl_matrix_get(H01, 2, risk1_index));
                 gsl_vector_add(TNN,TN);
                 gsl_vector_scale(TN, 1/gsl_matrix_get(H01, 2, risk1_index));
                 risk1_index--;
              }
              else
              {
                 for (j=j+1;j<k;j++)
                 {
                     for(t=0;t<p1a;t++)   gsl_matrix_set(D,t,t,gsl_matrix_get(FUNBSE,t,j));

                     if(p1a>1)
                     {
                         for(i=1;i<p1a;i++)
                         {
                             for(t=0;t<p1a-i;t++)   gsl_matrix_set(D,t,i+t,gsl_matrix_get(FUNBSE,p1a+t+(i-1)*(p1a-1),j));
                         }
                     }
                     for(t=0;t<p1a;t++)
                     {
                         for(i=0;i<t;i++)   gsl_matrix_set(D,t,i,gsl_matrix_get(D,i,t));
                     }

                     for (t=0;t<p1a;t++) gsl_vector_set(N,t,gsl_matrix_get(FUNBE,t,j));

                     for(u=0;u<p2;u++)
                     {
                         gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
                         gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));
                     }
                     gsl_matrix_scale(D,exp(MulVV(X,gammai)));
                     gsl_matrix_add(TD,D);
                     gsl_vector_scale(N,exp(MulVV(X,gammai)));
                     gsl_vector_add(TN,N);

                    if (j == k-1)
                    {
                       gsl_matrix_scale(TD, gsl_matrix_get(H01, 2, risk1_index));
                       gsl_matrix_add(TDD,TD);
                       gsl_matrix_scale(TD, 1/gsl_matrix_get(H01, 2, risk1_index));
                       gsl_vector_scale(TN, gsl_matrix_get(H01, 2, risk1_index));
                       gsl_vector_add(TNN,TN);
                       gsl_vector_scale(TN, 1/gsl_matrix_get(H01, 2, risk1_index));
                       risk1_index--;
                       break;
                    }
                    else if (gsl_matrix_get(C,j+1,0) != gsl_matrix_get(C,j,0))
                    {
                       gsl_matrix_scale(TD, gsl_matrix_get(H01, 2, risk1_index));
                       gsl_matrix_add(TDD,TD);
                       gsl_matrix_scale(TD, 1/gsl_matrix_get(H01, 2, risk1_index));
                       gsl_vector_scale(TN, gsl_matrix_get(H01, 2, risk1_index));
                       gsl_vector_add(TNN,TN);
                       gsl_vector_scale(TN, 1/gsl_matrix_get(H01, 2, risk1_index));
                       risk1_index--;
                       break;
                    }
                    else continue;
                 }
              }
            }
            else continue;
        }
        gsl_vector_set_zero(TN);
        for (j=0;j<k;j++)
        {
            if((int)gsl_matrix_get(C,j,1)==1)
            {
                for (t=0;t<p1a;t++) gsl_vector_set(N,t,gsl_matrix_get(FUNB,t,j));
                gsl_vector_add(TN, N);
            } else continue;
        }
        gsl_vector_sub(TN,TNN);

        status=inv_matrix(TDD);
        if(status==100)  return status;
        MulM(TDD,TN,N);


        for(j=0;j<p1a;j++)   gsl_vector_set(vee1,j,gsl_vector_get(vee1,j)+gsl_vector_get(N,j));

        gsl_matrix_set_zero(TD);
        gsl_matrix_set_zero(TDD);
        gsl_vector_set_zero(TN);
        gsl_vector_set_zero(TNN);

        gsl_vector_free(Z);
        gsl_vector_free(SZ);
        gsl_vector_free(X);
        gsl_vector_free(gammai);
        gsl_vector_free(SX);
        gsl_matrix_free(SZZ);
        gsl_matrix_free(ZZ);
        gsl_matrix_free(SXX);
        gsl_matrix_free(XX);

        gsl_matrix_free(FUNB);
        gsl_matrix_free(FUNE);
        gsl_matrix_free(FUNBSE);
        gsl_matrix_free(FUNBE);
        gsl_vector_free(N);
        gsl_vector_free(TN);
        gsl_matrix_free(D);
        gsl_matrix_free(TD);
        gsl_matrix_free(TDD);
        gsl_vector_free(TNN);
        gsl_vector_free(X_new);
        gsl_matrix_free(SXX_new);
        gsl_vector_free(SX_new);
        gsl_vector_free(SX_inter);
        gsl_vector_free(xtilde1);
        gsl_vector_free(xtilde);
        gsl_matrix_free(bs);
        gsl_vector_free(bi);
        return 0;
    }

    int GetCov(
               gsl_matrix *Cov,
               const gsl_vector *beta,
               const gsl_matrix *gamma,
               const gsl_vector *vee1,
               const gsl_matrix *H01,
               const double sigma,
               const gsl_matrix *sig,
               const gsl_matrix *Y,
               const gsl_matrix *C,
               const gsl_vector *M1,
               const gsl_matrix *Posbi,
               const gsl_matrix *Poscov,
               const int p1a,
               const int point,
               const std::vector<double> xs,
               const std::vector<double> ws
               )
    {

        int p1=beta->size;
        int p2=gamma->size2;
        int g =gamma->size1;
        int a =H01->size2;
        int d =Cov->size1;
        int k = M1->size;

        int i,p,q,j,t,u;

        double temp,temp1;

        gsl_vector *S = gsl_vector_calloc(d);
        gsl_matrix *SS= gsl_matrix_calloc(d,d);


        gsl_matrix *FUNB=gsl_matrix_calloc(p1a,k),
                   *FUNBS=gsl_matrix_calloc(p1a*(p1a+1)/2,k),
                   *FUNE=gsl_matrix_calloc(g,k),
                   *FUNBSE=gsl_matrix_calloc(g*p1a*(p1a+1)/2,k),
                   *FUNBE=gsl_matrix_calloc(g*p1a,k);

        int status;
        status = GetE(FUNB,FUNBS,FUNBSE,FUNBE,FUNE,beta,gamma,vee1,H01,sigma,sig,Y,C,M1,Posbi,Poscov,p1a,point,xs,ws);

        if (status==100) return status;

        gsl_matrix *VC = gsl_matrix_calloc(p1a,p1a);
        gsl_matrix *VI = gsl_matrix_calloc(p1a,p1a);
        gsl_matrix *HE = gsl_matrix_calloc(p1a,p1a);


        gsl_matrix_memcpy(VC,sig);
        if(p1a>1)
        {
            for(i=0;i<p1a;i++)
            {
                for(j=i+1;j<p1a;j++)    gsl_matrix_set(VC,j,i,gsl_matrix_get(VC,i,j));
            }
        }
        gsl_matrix_memcpy(VI,VC);
        gsl_permutation * vp = gsl_permutation_calloc(p1a);
        gsl_linalg_LU_decomp (VI, vp, &i);
        double vdet=gsl_linalg_LU_det(VI,i);

        printf("det=%f\n",vdet);

        //if(vdet<0.0005) return 100;

        status=inv_matrix(VC);
        if(status==100) return status;

        gsl_vector * Z = gsl_vector_calloc(p1);                         /* covariates for Y */
        gsl_vector * SZ= gsl_vector_calloc(p1);
        gsl_vector * xtilde = gsl_vector_calloc(p1a);
        gsl_vector * xtilde1 = gsl_vector_calloc(p1a);
        gsl_vector * bi = gsl_vector_calloc(p1a);
        gsl_matrix * bs = gsl_matrix_calloc(p1a,p1a);

        gsl_vector * X = gsl_vector_calloc(p2);                         /* covariates for C */
        gsl_vector * RX= gsl_vector_calloc(p2);
        /* covariates for subjects in risk set */
        gsl_vector *SX = gsl_vector_calloc(p2);
        gsl_vector *SX1 = gsl_vector_calloc(p2);
        gsl_vector *SRX= gsl_vector_calloc(p2);
        gsl_vector *SRXX= gsl_vector_calloc(p2);
        gsl_matrix *SXX1 = gsl_matrix_calloc(p2,a);
        gsl_matrix *SXX11 = gsl_matrix_calloc(p2,a);
        gsl_matrix *SRXX1 = gsl_matrix_calloc(p2,k);
        gsl_vector * gammai = gsl_vector_calloc(p2);

        gsl_vector *N=gsl_vector_calloc(p1a);
        gsl_vector *TN=gsl_vector_calloc(p1a);
        gsl_vector *TN1=gsl_vector_calloc(p1a);
        gsl_vector *TRN=gsl_vector_calloc(p1a);
        gsl_vector *TRNN=gsl_vector_calloc(p1a);
        gsl_matrix *TRNN1=gsl_matrix_calloc(p1a,k);
        gsl_matrix *TNN1=gsl_matrix_calloc(p1a,a);
        gsl_matrix *TNN11=gsl_matrix_calloc(p1a,a);



        p=0;
        gsl_matrix_set_zero(Cov);

        int risk1_index;
        int risk1_index_temp=a-1;
        int risk1_index_ttemp=a-1;
        int risk1_index_tttemp=a-1;
        int risk1_index_vtemp=a-1;
        int risk1_index_vttemp=a-1;
        int risk1_index_vtttemp=a-1;

        temp1=0;
        gsl_vector *CumuH01 = gsl_vector_calloc(a);
        for (j=0;j<a;j++)
        {
            temp1+=gsl_matrix_get(H01, 2, j);
            gsl_vector_set(CumuH01,j,temp1);
        }

        for(j=0;j<k;j++)
        {

            gsl_vector_set_zero(S);
            u=(int)gsl_vector_get(M1,j);

            /* calculate score for beta */

            gsl_vector_set_zero(SZ);
            for(q=p;q<u+p;q++)
            {
                for(t=0;t<p1a;t++)   gsl_vector_set(bi,t,gsl_matrix_get(FUNB,t,j));
                for(t=0;t<p1;t++)   gsl_vector_set(Z,t,gsl_matrix_get(Y,q,1+p1a+t));
                for(t=0;t<p1a;t++)   gsl_vector_set(xtilde,t,gsl_matrix_get(Y,q,1+t));
                temp = MulVV(beta,Z);

                gsl_vector_scale (Z,gsl_matrix_get(Y,q,0)-temp-MulVV(bi,xtilde));
                gsl_vector_add (SZ,Z);
            }

            gsl_vector_scale (SZ,1/sigma);
            for(q=0;q<p1;q++)  gsl_vector_set(S,q,gsl_vector_get(SZ,q));

            /* calculate score for sigma */

            temp=0;
            for(q=p;q<u+p;q++)
            {
                for(t=0;t<p1;t++)   gsl_vector_set(Z,t,gsl_matrix_get(Y,q,1+p1a+t));
                for(t=0;t<p1a;t++)   gsl_vector_set(xtilde,t,gsl_matrix_get(Y,q,1+t));
                for(t=0;t<p1a;t++)   gsl_matrix_set(bs,t,t,gsl_matrix_get(FUNBS,t,j));

                if(p1a>1)
                {
                    for(i=1;i<p1a;i++)
                    {
                        for(t=0;t<p1a-i;t++)   gsl_matrix_set(bs,t,i+t,gsl_matrix_get(FUNBS,p1a+t+(i-1)*(p1a-1),j));
                    }
                }

                for(t=0;t<p1a;t++)
                {
                    for(i=0;i<t;i++)   gsl_matrix_set(bs,t,i,gsl_matrix_get(bs,i,t));
                }

                MulM(bs,xtilde,xtilde1);

                temp+=gsl_pow_2(gsl_matrix_get(Y,q,0)-MulVV(Z,beta))-2*(gsl_matrix_get(Y,q,0)-MulVV(Z,beta))*MulVV(bi,xtilde)+MulVV(xtilde,xtilde1);
            }

            temp=temp/(2*sigma*sigma);
            temp=temp-(double)u/(2*sigma);
            gsl_vector_set(S,p1,temp);

            p=p+u;

            /* calculate score for gamma */
            /*  gamma11, gamma12 */
            for(u=0;u<p2;u++) gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));
            if (j == 0)
            {
                temp=0;
                risk1_index=risk1_index_temp;
                for (q=j;q<k;q++)
                {
                    for(u=0;u<p2;u++) gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));
                    temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q);
                    gsl_vector_scale(RX,exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q));
                    gsl_vector_add(SX,RX);

                    if (gsl_matrix_get(C,q,1) == 1)
                    {
                        if (q == k-1)
                        {
                            gsl_vector_scale(SX, gsl_matrix_get(H01,1,risk1_index)/(temp*temp));
                            gsl_matrix_set_col(SXX1,a-1-risk1_index,SX);
                            gsl_vector_add(SRXX,SX);
                            gsl_vector_scale(SX, 1/gsl_matrix_get(H01,1,risk1_index)*(temp*temp));
                            gsl_vector_scale(SX, 1/temp);
                            gsl_matrix_set_col(SXX11, a-1-risk1_index, SX);
                            gsl_vector_scale(SX, temp);
                            //gsl_matrix_set_col(SRXX1, j, SRX);
                            risk1_index--;
                        }
                        else if (gsl_matrix_get(C,q+1,0) != gsl_matrix_get(C,q,0))
                        {
                            gsl_vector_scale(SX, gsl_matrix_get(H01,1,risk1_index)/(temp*temp));
                            gsl_matrix_set_col(SXX1,a-1-risk1_index,SX);
                            gsl_vector_add(SRXX,SX);
                            gsl_vector_scale(SX, 1/gsl_matrix_get(H01,1,risk1_index)*(temp*temp));
                            gsl_vector_scale(SX, 1/temp);
                            gsl_matrix_set_col(SXX11, a-1-risk1_index, SX);
                            gsl_vector_scale(SX, temp);
                            //gsl_matrix_set_col(SRXX1, j, SRX);
                            risk1_index--;
                        }
                        else
                        {
                            for (q=q+1;q<k;q++)
                            {
                                for(u=0;u<p2;u++) gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));
                                temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q);
                                gsl_vector_scale(RX,exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q));
                                gsl_vector_add(SX,RX);
                                if (q == k-1)
                                {
                                    gsl_vector_scale(SX, gsl_matrix_get(H01,1,risk1_index)/(temp*temp));
                                    gsl_matrix_set_col(SXX1,a-1-risk1_index,SX);
                                    gsl_vector_add(SRXX,SX);
                                    gsl_vector_scale(SX, 1/gsl_matrix_get(H01,1,risk1_index)*(temp*temp));
                                    gsl_vector_scale(SX, 1/temp);
                                    gsl_matrix_set_col(SXX11, a-1-risk1_index, SX);
                                    gsl_vector_scale(SX, temp);
                                    //gsl_matrix_set_col(SRXX1, j, SRX);
                                    risk1_index--;
                                    break;
                                }
                                else if (gsl_matrix_get(C,q+1,0) != gsl_matrix_get(C,q,0))
                                {
                                    gsl_vector_scale(SX, gsl_matrix_get(H01,1,risk1_index)/(temp*temp));
                                    gsl_matrix_set_col(SXX1,a-1-risk1_index,SX);
                                    gsl_vector_add(SRXX,SX);
                                    gsl_vector_scale(SX, 1/gsl_matrix_get(H01,1,risk1_index)*(temp*temp));
                                    gsl_vector_scale(SX, 1/temp);
                                    gsl_matrix_set_col(SXX11, a-1-risk1_index, SX);
                                    gsl_vector_scale(SX, temp);
                                    //gsl_matrix_set_col(SRXX1, j, SRX);
                                    risk1_index--;
                                    break;
                                }
                                else continue;
                            }
                        }

                    }
                    else continue;
                }
                gsl_matrix_set_col(SRXX1,j,SRXX);
            }
            else
            {
                if (risk1_index_temp>=0)
                {
                    if (gsl_matrix_get(C,j,0) >= gsl_matrix_get(H01,0,risk1_index_temp))
                    {
                        gsl_matrix_get_col(SRXX,SRXX1,j-1);
                        gsl_matrix_set_col(SRXX1,j,SRXX);
                    }
                    else
                    {
                        risk1_index_temp--;
                        if (risk1_index_temp>=0)
                        {
                            gsl_matrix_get_col(SX1, SXX1, a-1-risk1_index_temp-1);
                            gsl_matrix_get_col(SRXX,SRXX1,j-1);
                            gsl_vector_sub(SRXX, SX1);
                            gsl_matrix_set_col(SRXX1,j,SRXX);
                        }
                    }
                }
                else
                {
                    gsl_matrix_get_col(SRX, SRXX1, j);
                    risk1_index_temp=0;
                }
            }

            gsl_matrix_get_col(SRX, SRXX1, j);

            if (j==0)
            {
                for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
                gsl_vector_scale(SX, gsl_vector_get(CumuH01, risk1_index_ttemp));
                //gsl_matrix_get_col(SRX, SRXX1, j);
                gsl_vector_sub(SRX,SX);
            }
            else if (gsl_matrix_get(C,j,0) >= gsl_matrix_get(H01,0,risk1_index_ttemp))
            {
                //gsl_matrix_get_col(SRX,SRXX1,j);
                for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
                gsl_vector_scale(SX, gsl_vector_get(CumuH01, risk1_index_ttemp));
                gsl_vector_sub(SRX,SX);
            }
            else
            {
                risk1_index_ttemp--;
                if (risk1_index_ttemp>=0)
                {
                    //gsl_matrix_get_col(SRX,SRXX1,j);
                    for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
                    gsl_vector_scale(SX, gsl_vector_get(CumuH01, risk1_index_ttemp));
                    gsl_vector_sub(SRX,SX);
                }
                else
                {
                    gsl_matrix_get_col(SRX,SRXX1,j);
                    risk1_index_ttemp=0;
                }
            }

            for(u=0;u<p2;u++) gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));

            gsl_vector_scale(SRX, exp(MulVV(SX,gammai))*gsl_matrix_get(FUNE,0,j));

            if (gsl_matrix_get(C,j,0) >= gsl_matrix_get(H01,0,risk1_index_tttemp))
            {
                if ((int)gsl_matrix_get(C,j,1) == 1)
                {
                    for(u=0;u<p2;u++)  gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
                    gsl_matrix_get_col(SX,SXX11,a-1-risk1_index_tttemp);
                    gsl_vector_sub(X,SX);
                    gsl_vector_add(X,SRX);



                    for(q=0;q<p2;q++)   gsl_vector_set(S,q+p1+1,gsl_vector_get(X,q));



                }
                else
                {
                    for(q=0;q<p2;q++)   gsl_vector_set(S,q+p1+1,gsl_vector_get(SRX,q));
                }
            }
            else
            {
                risk1_index_tttemp--;
                if (risk1_index_tttemp>=0)
                {
                    if ((int)gsl_matrix_get(C,j,1) == 1)
                    {
                        for(u=0;u<p2;u++)  gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
                        gsl_matrix_get_col(SX,SXX11,a-1-risk1_index_tttemp);
                        gsl_vector_sub(X,SX);
                        gsl_vector_add(X,SRX);
                        for(q=0;q<p2;q++)   gsl_vector_set(S,q+p1+1,gsl_vector_get(X,q));

                    }
                    else
                    {
                        for(q=0;q<p2;q++)   gsl_vector_set(S,q+p1+1,gsl_vector_get(SRX,q));
                    }
                }
                else
                {
                    risk1_index_tttemp=0;
                    for(q=0;q<p2;q++)   gsl_vector_set(S,q+p1+1,gsl_vector_get(SRX,q));
                }
            }

            /* calculate score for vee */
            /*  vee1 */
            for(u=0;u<p2;u++) gsl_vector_set(gammai,u,gsl_matrix_get(gamma,0,u));
            if (j == 0)
            {
                temp=0;
                gsl_vector_set_zero(TN);
                gsl_vector_set_zero(TRN);
                risk1_index=risk1_index_vtemp;
                for (q=j;q<k;q++)
                {
                    for(t=0;t<p1a;t++) gsl_vector_set(N,t,gsl_matrix_get(FUNBE,t,q));
                    for(u=0;u<p2;u++) gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));
                    temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q);
                    gsl_vector_scale(N,exp(MulVV(RX,gammai)));
                    gsl_vector_add(TN,N);

                    if (gsl_matrix_get(C,q,1) == 1)
                    {
                        if (q == k-1)
                        {
                            gsl_vector_scale(TN, gsl_matrix_get(H01,1,risk1_index)/(temp*temp));
                            gsl_matrix_set_col(TNN1,a-1-risk1_index,TN);
                            gsl_vector_add(TRNN,TN);
                            gsl_vector_scale(TN, 1/gsl_matrix_get(H01,1,risk1_index)*(temp*temp));
                            gsl_vector_scale(TN, 1/temp);
                            gsl_matrix_set_col(TNN11, a-1-risk1_index, TN);
                            gsl_vector_scale(TN, temp);
                            risk1_index--;
                        }
                        else if (gsl_matrix_get(C,q+1,0) != gsl_matrix_get(C,q,0))
                        {
                            gsl_vector_scale(TN, gsl_matrix_get(H01,1,risk1_index)/(temp*temp));
                            gsl_matrix_set_col(TNN1,a-1-risk1_index,TN);
                            gsl_vector_add(TRNN,TN);
                            gsl_vector_scale(TN, 1/gsl_matrix_get(H01,1,risk1_index)*(temp*temp));
                            gsl_vector_scale(TN, 1/temp);
                            gsl_matrix_set_col(TNN11, a-1-risk1_index, TN);
                            gsl_vector_scale(TN, temp);
                            risk1_index--;
                        }
                        else
                        {
                            for (q=q+1;q<k;q++)
                            {
                                for(t=0;t<p1a;t++) gsl_vector_set(N,t,gsl_matrix_get(FUNBE,t,q));
                                for(u=0;u<p2;u++) gsl_vector_set(RX,u,gsl_matrix_get(C,q,2+u));
                                temp+=exp(MulVV(RX,gammai))*gsl_matrix_get(FUNE,0,q);
                                gsl_vector_scale(N,exp(MulVV(RX,gammai)));
                                gsl_vector_add(TN,N);
                                if (q == k-1)
                                {
                                    gsl_vector_scale(TN, gsl_matrix_get(H01,1,risk1_index)/(temp*temp));
                                    gsl_matrix_set_col(TNN1,a-1-risk1_index,TN);
                                    gsl_vector_add(TRNN,TN);
                                    gsl_vector_scale(TN, 1/gsl_matrix_get(H01,1,risk1_index)*(temp*temp));
                                    gsl_vector_scale(TN, 1/temp);
                                    gsl_matrix_set_col(TNN11, a-1-risk1_index, TN);
                                    gsl_vector_scale(TN, temp);
                                    risk1_index--;
                                    break;
                                }
                                else if (gsl_matrix_get(C,q+1,0) != gsl_matrix_get(C,q,0))
                                {
                                    gsl_vector_scale(TN, gsl_matrix_get(H01,1,risk1_index)/(temp*temp));
                                    gsl_matrix_set_col(TNN1,a-1-risk1_index,TN);
                                    gsl_vector_add(TRNN,TN);
                                    gsl_vector_scale(TN, 1/gsl_matrix_get(H01,1,risk1_index)*(temp*temp));
                                    gsl_vector_scale(TN, 1/temp);
                                    gsl_matrix_set_col(TNN11, a-1-risk1_index, TN);
                                    gsl_vector_scale(TN, temp);
                                    risk1_index--;
                                    break;
                                }
                                else continue;
                            }
                        }

                    }
                    else continue;
                }
                gsl_matrix_set_col(TRNN1,j,TRNN);
            }
            else
            {
                if (risk1_index_vtemp>=0)
                {
                    if (gsl_matrix_get(C,j,0) >= gsl_matrix_get(H01,0,risk1_index_vtemp))
                    {
                        gsl_matrix_get_col(TRNN,TRNN1,j-1);
                        gsl_matrix_set_col(TRNN1,j,TRNN);
                    }
                    else
                    {
                        risk1_index_vtemp--;
                        if (risk1_index_vtemp>=0)
                        {
                            gsl_matrix_get_col(TN1, TNN1, a-1-risk1_index_vtemp-1);
                            gsl_matrix_get_col(TRNN,TRNN1,j-1);
                            gsl_vector_sub(TRNN, TN1);
                            gsl_matrix_set_col(TRNN1,j,TRNN);
                        }
                    }
                }
                else
                {
                    gsl_matrix_get_col(TRN, TRNN1, j);
                    risk1_index_vtemp=0;
                }
            }

            gsl_matrix_get_col(TRN, TRNN1, j);
            for(u=0;u<p2;u++)  gsl_vector_set(SX,u,gsl_matrix_get(C,j,2+u));
            gsl_vector_scale(TRN,gsl_matrix_get(FUNE,0,j)*exp(MulVV(SX,gammai)));

            if (j==0)
            {
                for (t=0;t<p1a;t++) gsl_vector_set(N,t,gsl_matrix_get(FUNBE,t,j));
                gsl_vector_scale(N, gsl_vector_get(CumuH01, risk1_index_vttemp)*exp(MulVV(SX,gammai)));
                gsl_vector_sub(TRN,N);
            }
            else if (gsl_matrix_get(C,j,0) >= gsl_matrix_get(H01,0,risk1_index_vttemp))
            {
                for (t=0;t<p1a;t++) gsl_vector_set(N,t,gsl_matrix_get(FUNBE,t,j));
                gsl_vector_scale(N, gsl_vector_get(CumuH01, risk1_index_vttemp)*exp(MulVV(SX,gammai)));
                gsl_vector_sub(TRN,N);
            }
            else
            {
                risk1_index_vttemp--;
                if (risk1_index_vttemp>=0)
                {
                    for (t=0;t<p1a;t++) gsl_vector_set(N,t,gsl_matrix_get(FUNBE,t,j));
                    gsl_vector_scale(N, gsl_vector_get(CumuH01, risk1_index_vttemp)*exp(MulVV(SX,gammai)));
                    gsl_vector_sub(TRN,N);
                }
                else
                {
                    risk1_index_vttemp=0;
                }
            }


            if (gsl_matrix_get(C,j,0) >= gsl_matrix_get(H01,0,risk1_index_vtttemp))
            {
                if ((int)gsl_matrix_get(C,j,1) == 1)
                {
                    for (t=0;t<p1a;t++) gsl_vector_set(N,t,gsl_matrix_get(FUNBE,t,j));
                    for(u=0;u<p2;u++)  gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
                    gsl_matrix_get_col(TN,TNN11,a-1-risk1_index_vtttemp);
                    for (t=0;t<p1a;t++) gsl_vector_set(TN,t,gsl_matrix_get(FUNB,t,j)-gsl_vector_get(TN,t));
                    gsl_vector_add(TN,TRN);
                    for(u=0;u<p1a;u++)  gsl_vector_set(S,p1+p2+1+u,gsl_vector_get(TN,u));
                }
                else
                {
                    for(u=0;u<p1a;u++)  gsl_vector_set(S,p1+p2+1+u,gsl_vector_get(TRN,u));
                }
            }
            else
            {
                risk1_index_vtttemp--;
                if (risk1_index_vtttemp>=0)
                {
                    if ((int)gsl_matrix_get(C,j,1) == 1)
                    {
                        for (t=0;t<p1a;t++) gsl_vector_set(N,t,gsl_matrix_get(FUNBE,t,j));
                        for(u=0;u<p2;u++)  gsl_vector_set(X,u,gsl_matrix_get(C,j,2+u));
                        gsl_matrix_get_col(TN,TNN11,a-1-risk1_index_vtttemp);
                        for (t=0;t<p1a;t++) gsl_vector_set(TN,t,gsl_matrix_get(FUNB,t,j)-gsl_vector_get(TN,t));
                        gsl_vector_add(TN,TRN);
                        for(u=0;u<p1a;u++)  gsl_vector_set(S,p1+p2+1+u,gsl_vector_get(TN,u));
                    }
                    else
                    {
                        for(u=0;u<p1a;u++)  gsl_vector_set(S,p1+p2+1+u,gsl_vector_get(TRN,u));
                    }
                }
                else
                {
                    risk1_index_vtttemp=0;
                    for(u=0;u<p1a;u++)  gsl_vector_set(S,p1+p2+1+u,gsl_vector_get(TRN,u));
                }
            }

            /*  Sigma matrix  */


            for(t=0;t<p1a;t++)   gsl_matrix_set(bs,t,t,gsl_matrix_get(FUNBS,t,j));
            if(p1a>1)
            {
                for(i=1;i<p1a;i++)
                {
                    for(t=0;t<p1a-i;t++)   gsl_matrix_set(bs,t,i+t,gsl_matrix_get(FUNBS,p1a+t+(i-1)*(p1a-1),j));
                }

                for(t=0;t<p1a;t++)
                {
                    for(i=0;i<t;i++)   gsl_matrix_set(bs,t,i,gsl_matrix_get(bs,i,t));
                }
            }


            gsl_matrix_memcpy(VI,bs);

            MulMM(VC,VI,HE);
            MulMM(HE,VC,VI);

            if (p1a==1)
            {
                gsl_vector_set(S,p1+1+p2+p1a,(gsl_matrix_get(VI,0,0)-gsl_matrix_get(VC,0,0))/2);
            } else
            {
                for (t=0;t<p1a;t++)
                {
                    gsl_vector_set(S,p1+1+p2+p1a+t,(gsl_matrix_get(VI,t,t)-gsl_matrix_get(VC,t,t))/2);
                }
                for(q=1;q<p1a;q++)
                {
                    for(t=0;t<p1a-q;t++) gsl_vector_set(S,p1+1+p2+p1a+p1a+t+(q-1)*(p1a-1),(gsl_matrix_get(VI,t,q+t)-gsl_matrix_get(VC,t,q+t)));
                }
            }
            MulV(S,SS);
            gsl_matrix_add(Cov,SS);
        }


        gsl_vector_free(Z);
        gsl_vector_free(SZ);
        gsl_vector_free(X);
        gsl_vector_free(RX);
        gsl_vector_free(SX);
        gsl_vector_free(SRX);
        gsl_vector_free(gammai);
        gsl_vector_free(S);
        gsl_matrix_free(SS);

        gsl_vector_free(xtilde);
        gsl_vector_free(xtilde1);
        gsl_vector_free(bi);
        gsl_matrix_free(bs);

        gsl_vector_free(SX1);
        gsl_vector_free(SRXX);
        gsl_matrix_free(SXX1);
        gsl_matrix_free(SXX11);
        gsl_matrix_free(SRXX1);
        gsl_matrix_free(FUNB);
        gsl_matrix_free(FUNBS);
        gsl_matrix_free(FUNE);
        gsl_matrix_free(FUNBSE);
        gsl_matrix_free(FUNBE);

        gsl_vector_free(N);
        gsl_vector_free(TN);
        gsl_vector_free(TRN);

        gsl_vector_free(TN1);
        gsl_vector_free(TRNN);

        gsl_matrix_free(TRNN1);
        gsl_matrix_free(TNN1);
        gsl_matrix_free(TNN11);

        gsl_matrix_free(VC);
        gsl_matrix_free(VI);
        gsl_matrix_free(HE);
        gsl_permutation_free(vp);


        return 0;
    }

    int GetE(
             gsl_matrix *FUNB,
             gsl_matrix *FUNBS,
             gsl_matrix *FUNBSE,
             gsl_matrix *FUNBE,
             gsl_matrix *FUNE,
             const gsl_vector *beta,
             const gsl_matrix *gamma,
             const gsl_vector *vee1,
             const gsl_matrix *H01,
             const double sigma,
             const gsl_matrix *sig,
             const gsl_matrix *Y,
             const gsl_matrix *C,
             const gsl_vector *M1,
             const gsl_matrix *Posbi,
             const gsl_matrix *Poscov,
             const int p1a,
             const int point,
             const std::vector<double> xs,
             const std::vector<double> ws
             )
    {

        gsl_matrix_set_zero(FUNB);
        gsl_matrix_set_zero(FUNBS);
        gsl_matrix_set_zero(FUNBSE);
        gsl_matrix_set_zero(FUNBE);

        gsl_matrix_set_zero(FUNE);


        int p1=beta->size;
        int p2=gamma->size2;
        int a =H01->size2;

        int k = M1->size;


        int i,j,q,t,m;
        double mu,dem,temp;
        double cuh01,haz01,xgamma1;


        gsl_vector *Z = gsl_vector_calloc(p1),
                   *X = gsl_vector_calloc(p2),
                   *xtilde = gsl_vector_calloc(p1a),
                   *gammai = gsl_vector_calloc(p2);

        int db0,db1,db2;

        gsl_vector *xi = gsl_vector_calloc(point);
        gsl_vector *wi = gsl_vector_calloc(point);


        for(i=0;i<point/2;i++)   gsl_vector_set(xi,i,xs[i]);
        for(i=0;i<point/2;i++)   gsl_vector_set(wi,i,ws[i]);

        for(i=0;i<point/2;i++)   gsl_vector_set(xi,point-i-1, 0-gsl_vector_get(xi,i));
        for(i=0;i<point/2;i++)   gsl_vector_set(wi,point-i-1, gsl_vector_get(wi,i));

        gsl_vector *ti=gsl_vector_calloc(p1a);
        gsl_vector *bi=gsl_vector_calloc(p1a);
        gsl_vector *ci=gsl_vector_calloc(p1a);
        gsl_matrix *covi=gsl_matrix_calloc(p1a,p1a);
        gsl_vector *rii=gsl_vector_calloc(p1a);
        gsl_vector *tii=gsl_vector_calloc(p1a);

        m=0;

        gsl_vector *CUH01 = gsl_vector_calloc(k);
        gsl_vector *HAZ01 = gsl_vector_calloc(k);
        int risk1_index = a-1;

        double temp1=0;
        gsl_vector *CumuH01 = gsl_vector_calloc(a);
        for (j=0;j<a;j++)
        {
            temp1+=gsl_matrix_get(H01, 2, j);
            gsl_vector_set(CumuH01,j,temp1);
        }

        for (j=0;j<k;j++)
        {
            if (risk1_index>=0)
            {
                if (gsl_matrix_get(C,j,0) >= gsl_matrix_get(H01,0,risk1_index))
                {
                    gsl_vector_set(CUH01,j,gsl_vector_get(CumuH01,risk1_index));
                }
                else
                {
                    risk1_index--;
                    if (risk1_index>=0)
                    {
                        gsl_vector_set(CUH01,j,gsl_vector_get(CumuH01,risk1_index));
                    }
                }
            }
            else
            {
                risk1_index=0;
            }
        }

        risk1_index = a-1;

        for (j=0;j<k;j++)
        {
            if (risk1_index>=0)
            {
                   //gsl_vector_set(CUH01,j,gsl_vector_get(CumuH01,risk1_index));
                   if (gsl_matrix_get(C,j,0) == gsl_matrix_get(H01,0,risk1_index))
                   {
                       gsl_vector_set(HAZ01,j,gsl_matrix_get(H01,2,risk1_index));
                   }
                   if ((int)gsl_matrix_get(C,j,1) == 1)
                   {
                       if (j == k-1)
                       {
                           risk1_index--;
                       }
                       else if (gsl_matrix_get(C,j+1,0) != gsl_matrix_get(C,j,0))
                       {
                           risk1_index--;
                       }
                       else
                       {
                           for (j=j+1;j<k;j++)
                           {
                               //gsl_vector_set(CUH01,j,gsl_vector_get(CumuH01,risk1_index));
                               if (gsl_matrix_get(C,j,0) == gsl_matrix_get(H01,0,risk1_index))
                               {
                                   gsl_vector_set(HAZ01,j,gsl_matrix_get(H01,2,risk1_index));
                               }
                               if (j == k-1)
                               {
                                   risk1_index--;
                                   break;
                               }
                               else if (gsl_matrix_get(C,j+1,0) != gsl_matrix_get(C,j,0))
                               {
                                   risk1_index--;
                                   break;
                               }
                               else continue;
                           }
                       }
                   }
                   else continue;
            }
            else continue;
        }

        gsl_matrix *VC = gsl_matrix_calloc(p1a,p1a);
        gsl_matrix_memcpy(VC,sig);
        if(p1a>1)
        {
            for(i=0;i<p1a;i++)
            {
                for(j=i+1;j<p1a;j++)    gsl_matrix_set(VC,j,i,gsl_matrix_get(VC,i,j));
            }
        }

        int status;
        status=inv_matrix(VC);
        gsl_vector *tiii = gsl_vector_calloc(p1a);
        double u;

        for(j=0;j<k;j++)
        {
            dem=0;
            q=(int)gsl_vector_get(M1,j);
            cuh01=gsl_vector_get(CUH01,j);
            haz01=gsl_vector_get(HAZ01,j);
            for(i=0;i<p2;i++)
            {
                gsl_vector_set(X,i,gsl_matrix_get(C,j,2+i));
                gsl_vector_set(gammai,i,gsl_matrix_get(gamma,0,i));
            }
            xgamma1=MulVV(X,gammai);

            /*Extract Bayes estimate*/
            for (i=0;i<p1a;i++) gsl_vector_set(bi, i, gsl_matrix_get(Posbi, j, i));
            for (i=0;i<p1a;i++)
            {
                gsl_matrix_get_row(ci, Poscov, j*p1a+i);
                gsl_matrix_set_row(covi, i, ci);
            }


            if(p1a==1)
            {

                for(db0=0;db0<point;db0++)
                {
                        gsl_vector_set(rii,0,gsl_vector_get(xi, db0));
                        MulM(covi, rii, tii);
                        gsl_vector_scale(tii, sqrt(2));
                        gsl_vector_add(tii, bi);
                        gsl_vector_memcpy(ti, tii);

                        temp=1;

                        for(i=0;i<q;i++)
                        {
                            for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                            for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));
                            mu=MulVV(Z,beta);
                        temp*=exp(-1/(2*sigma)*gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-MulVV(ti,xtilde)));
                        }

                        if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+MulVV(vee1,ti));

                        temp*=exp(0-cuh01*exp(xgamma1+MulVV(vee1,ti)));
                        temp*=gsl_vector_get(wi,db0);
                        gsl_blas_dgemv(CblasNoTrans, 1.0, VC, tii, 0.0, tiii);
                        gsl_blas_ddot(tii, tiii, &u);
                        temp*=exp(gsl_pow_2(gsl_blas_dnrm2(rii)) - u/2);

                        dem+=temp;
                        for (i=0;i<p1a;i++)
                        {
                           gsl_matrix_set(FUNB,i,j,gsl_matrix_get(FUNB,i,j)+temp*gsl_vector_get(ti,i));
                           gsl_matrix_set(FUNBS,i,j,gsl_matrix_get(FUNBS,i,j)+temp*gsl_pow_2(gsl_vector_get(ti,i)));
                        }

                        gsl_matrix_set(FUNE,0,j,gsl_matrix_get(FUNE,0,j)+temp*exp(MulVV(vee1,ti)));

                        for (i=0;i<p1a;i++)
                        {
                            gsl_matrix_set(FUNBSE,i,j,gsl_matrix_get(FUNBSE,i,j)+temp*gsl_pow_2(gsl_vector_get(ti,i))*exp(MulVV(vee1,ti)));
                        }

                        for (i=0;i<p1a;i++)
                        {
                            gsl_matrix_set(FUNBE, i, j, gsl_matrix_get(FUNBE,i,j)+temp*gsl_vector_get(ti,i)*exp(MulVV(vee1,ti)));
                        }

                }
            }

            if(p1a==2)
            {
                for (db0=0;db0<point;db0++)
                {
                    for(db1=0;db1<point;db1++)
                    {
                            gsl_vector_set(rii,0,gsl_vector_get(xi, db0));
                            gsl_vector_set(rii,1,gsl_vector_get(xi, db1));
                            MulM(covi, rii, tii);
                            gsl_vector_scale(tii, sqrt(2));
                            gsl_vector_add(tii, bi);
                            gsl_vector_memcpy(ti, tii);
                            temp=1;

                            for(i=0;i<q;i++)
                            {
                                for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                                for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));
                                mu=MulVV(Z,beta);
                                temp*=exp(-1/(2*sigma)*gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-MulVV(ti,xtilde)));
                            }

                            if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+MulVV(vee1,ti));

                            temp*=exp(0-cuh01*exp(xgamma1+MulVV(vee1,ti)));
                            temp*=gsl_vector_get(wi,db0)*gsl_vector_get(wi,db1);
                            gsl_blas_dgemv(CblasNoTrans, 1.0, VC, tii, 0.0, tiii);
                            gsl_blas_ddot(tii, tiii, &u);
                            temp*=exp(gsl_pow_2(gsl_blas_dnrm2(rii)) - u/2);
                            dem+=temp;
                            for (i=0;i<p1a;i++)
                            {
                               gsl_matrix_set(FUNB,i,j,gsl_matrix_get(FUNB,i,j)+temp*gsl_vector_get(ti,i));
                               gsl_matrix_set(FUNBS,i,j,gsl_matrix_get(FUNBS,i,j)+temp*gsl_pow_2(gsl_vector_get(ti,i)));
                            }


                            for(i=1;i<p1a;i++)
                            {
                                for(t=0;t<p1a-i;t++)   gsl_matrix_set(FUNBS,p1a+t+(i-1)*(p1a-1),j,gsl_matrix_get(FUNBS,p1a+t+(i-1)*(p1a-1),j)
                                                      +temp*gsl_vector_get(ti,t)*gsl_vector_get(ti,t+i));
                            }

                            gsl_matrix_set(FUNE,0,j,gsl_matrix_get(FUNE,0,j)+temp*exp(MulVV(vee1,ti)));

                            for (i=0;i<p1a;i++)
                            {
                                gsl_matrix_set(FUNBSE,i,j,gsl_matrix_get(FUNBSE,i,j)+temp*gsl_pow_2(gsl_vector_get(ti,i))*exp(MulVV(vee1,ti)));
                            }

                            for(i=1;i<p1a;i++)
                            {
                                for(t=0;t<p1a-i;t++)
                                {
                                    gsl_matrix_set(FUNBSE,p1a+t+(i-1)*(p1a-1),j,gsl_matrix_get(FUNBSE,p1a+t+(i-1)*(p1a-1),j)
                                                          +temp*gsl_vector_get(ti,t)*gsl_vector_get(ti,t+i)*exp(MulVV(vee1,ti)));
                                }
                            }

                            for (i=0;i<p1a;i++)
                            {
                                gsl_matrix_set(FUNBE, i, j, gsl_matrix_get(FUNBE,i,j)+temp*gsl_vector_get(ti,i)*exp(MulVV(vee1,ti)));
                            }

                    }
                }

            }



            if(p1a==3)
            {
                for (db0=0;db0<point;db0++)
                {
                    for (db1=0;db1<point;db1++)
                    {
                        for(db2=0;db2<point;db2++)
                        {

                                gsl_vector_set(rii,0,gsl_vector_get(xi, db0));
                                gsl_vector_set(rii,1,gsl_vector_get(xi, db1));
                                gsl_vector_set(rii,2,gsl_vector_get(xi, db2));
                                MulM(covi, rii, tii);
                                gsl_vector_scale(tii, sqrt(2));
                                gsl_vector_add(tii, bi);
                                gsl_vector_memcpy(ti, tii);

                                temp=1;

                                for(i=0;i<q;i++)
                                {
                                    for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                                    for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));
                                    mu=MulVV(Z,beta);
                                temp*=exp(-1/(2*sigma)*gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-MulVV(ti,xtilde)));
                                }

                                if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+MulVV(vee1,ti));

                                temp*=exp(0-cuh01*exp(xgamma1+MulVV(vee1,ti)));
                                temp*=gsl_vector_get(wi,db0)*gsl_vector_get(wi,db1)*gsl_vector_get(wi,db2);
                                gsl_blas_dgemv(CblasNoTrans, 1.0, VC, tii, 0.0, tiii);
                                gsl_blas_ddot(tii, tiii, &u);
                                temp*=exp(gsl_pow_2(gsl_blas_dnrm2(rii)) - u/2);

                                dem+=temp;
                                for (i=0;i<p1a;i++)
                                {
                                   gsl_matrix_set(FUNB,i,j,gsl_matrix_get(FUNB,i,j)+temp*gsl_vector_get(ti,i));
                                   gsl_matrix_set(FUNBS,i,j,gsl_matrix_get(FUNBS,i,j)+temp*gsl_pow_2(gsl_vector_get(ti,i)));
                                }

                                for(i=1;i<p1a;i++)
                                {
                                    for(t=0;t<p1a-i;t++)   gsl_matrix_set(FUNBS,p1a+t+(i-1)*(p1a-1),j,gsl_matrix_get(FUNBS,p1a+t+(i-1)*(p1a-1),j)
                                                          +temp*gsl_vector_get(ti,t)*gsl_vector_get(ti,t+i));
                                }

                                gsl_matrix_set(FUNE,0,j,gsl_matrix_get(FUNE,0,j)+temp*exp(MulVV(vee1,ti)));

                                for (i=0;i<p1a;i++)
                                {
                                    gsl_matrix_set(FUNBSE,i,j,gsl_matrix_get(FUNBSE,i,j)+temp*gsl_pow_2(gsl_vector_get(ti,i))*exp(MulVV(vee1,ti)));
                                }

                                for(i=1;i<p1a;i++)
                                {
                                    for(t=0;t<p1a-i;t++)
                                    {
                                        gsl_matrix_set(FUNBSE,p1a+t+(i-1)*(p1a-1),j,gsl_matrix_get(FUNBSE,p1a+t+(i-1)*(p1a-1),j)
                                                              +temp*gsl_vector_get(ti,t)*gsl_vector_get(ti,t+i)*exp(MulVV(vee1,ti)));
                                    }
                                }

                                for (i=0;i<p1a;i++)
                                {
                                    gsl_matrix_set(FUNBE, i, j, gsl_matrix_get(FUNBE,i,j)+temp*gsl_vector_get(ti,i)*exp(MulVV(vee1,ti)));
                                }

                        }
                    }
                }


            }

            if(dem==0) return 100;

            for(i=0;i<p1a;i++) gsl_matrix_set(FUNB,i,j,gsl_matrix_get(FUNB,i,j)/dem);

            for(i=0;i<p1a*(p1a+1)/2;i++) gsl_matrix_set(FUNBS,i,j,gsl_matrix_get(FUNBS,i,j)/dem);

            gsl_matrix_set(FUNE,0,j,gsl_matrix_get(FUNE,0,j)/dem);

            for (i=0;i<p1a*(p1a+1)/2;i++) gsl_matrix_set(FUNBSE,i,j,gsl_matrix_get(FUNBSE,i,j)/dem);

            for(i=0;i<p1a;i++) gsl_matrix_set(FUNBE,i,j,gsl_matrix_get(FUNBE,i,j)/dem);

            m+=q;

        }


        gsl_vector_free(Z);
        gsl_vector_free(X);
        gsl_vector_free(xtilde);
        gsl_vector_free(gammai);


        gsl_vector_free(xi);
        gsl_vector_free(ti);
        gsl_vector_free(wi);
        gsl_vector_free(bi);
        gsl_matrix_free(covi);
        gsl_matrix_free(VC);
        gsl_vector_free(ci);
        gsl_vector_free(rii);
        gsl_vector_free(tii);
        gsl_vector_free(tiii);
        gsl_vector_free(CUH01);
        gsl_vector_free(HAZ01);
        gsl_vector_free(CumuH01);



        return 0;
    }

    double Getloglik(
                     const gsl_vector *beta,
                     const gsl_matrix *gamma,
                     const gsl_vector *vee1,
                     const gsl_matrix *H01,
                     const double sigma,
                     const gsl_matrix *sig,
                     const gsl_matrix *Y,
                     const gsl_matrix *C,
                     const gsl_vector *M1,
                     const gsl_matrix *Posbi,
                     const gsl_matrix *Poscov,
                     const int p1a,
                     const int point,
                     const std::vector<double> xs,
                     const std::vector<double> ws
                     )

    {
        int p1=beta->size;
        int p2=gamma->size2;
        int a =H01->size2;
        int k = M1->size;

        int i,j,q,t,m;
        double mu,temp1,temp;
        double cuh01,haz01,xgamma1;


        gsl_vector *Z = gsl_vector_calloc(p1),
                   *X = gsl_vector_calloc(p2),
                   *xtilde = gsl_vector_calloc(p1a),
                   *gammai = gsl_vector_calloc(p2);



        int db0,db1,db2;


        gsl_vector *xi = gsl_vector_calloc(point);
        gsl_vector *wi = gsl_vector_calloc(point);


        for(i=0;i<point/2;i++)   gsl_vector_set(xi,i,xs[i]);
        for(i=0;i<point/2;i++)   gsl_vector_set(wi,i,ws[i]);

        for(i=0;i<point/2;i++)   gsl_vector_set(xi,point-i-1, 0-gsl_vector_get(xi,i));
        for(i=0;i<point/2;i++)   gsl_vector_set(wi,point-i-1, gsl_vector_get(wi,i));


        gsl_vector *ti=gsl_vector_calloc(p1a);
        gsl_vector *bi=gsl_vector_calloc(p1a);
        gsl_vector *ci=gsl_vector_calloc(p1a);
        gsl_matrix *covi=gsl_matrix_calloc(p1a,p1a);
        gsl_vector *rii=gsl_vector_calloc(p1a);
        gsl_vector *tii=gsl_vector_calloc(p1a);

        gsl_matrix *VC = gsl_matrix_calloc(p1a,p1a);
        gsl_matrix_memcpy(VC,sig);
        if(p1a>1)
        {
            for(i=0;i<p1a;i++)
            {
                for(j=i+1;j<p1a;j++)    gsl_matrix_set(VC,j,i,gsl_matrix_get(VC,i,j));
            }
        }

        int status;
        status=inv_matrix(VC);
        gsl_vector *tiii = gsl_vector_calloc(p1a);
        double uu;

        double loglik=0;

        int risk1_index = a-1;
        gsl_vector *CumuH01 = gsl_vector_calloc(a);
        gsl_vector *CUH01 = gsl_vector_calloc(k);
        gsl_vector *HAZ01 = gsl_vector_calloc(k);
        temp1=0;
        for (j=0;j<a;j++)
        {
            temp1+=gsl_matrix_get(H01, 2, j);
            gsl_vector_set(CumuH01,j,temp1);
        }

        for (j=0;j<k;j++)
        {
            if (risk1_index>=0)
            {
                   gsl_vector_set(CUH01,j,gsl_vector_get(CumuH01,risk1_index));
                   if (gsl_matrix_get(C,j,1) == 1)
                   {

                       if (j == k-1)
                       {
                           gsl_vector_set(HAZ01,j,gsl_matrix_get(H01,2,risk1_index));
                           risk1_index--;
                       }
                       else if (gsl_matrix_get(C,j+1,0) != gsl_matrix_get(C,j,0))
                       {
                           gsl_vector_set(HAZ01,j,gsl_matrix_get(H01,2,risk1_index));
                           risk1_index--;
                       }

                       else
                       {
                           for (j=j+1;j<k;j++)
                           {
                               gsl_vector_set(CUH01,j,gsl_vector_get(CumuH01,risk1_index));
                               if (j == k-1)
                               {
                                   gsl_vector_set(HAZ01,j,gsl_matrix_get(H01,2,risk1_index));
                                   risk1_index--;
                                   break;
                               }
                               else if (gsl_matrix_get(C,j+1,0) != gsl_matrix_get(C,j,0))
                               {
                                   gsl_vector_set(HAZ01,j,gsl_matrix_get(H01,2,risk1_index));
                                   risk1_index--;
                                   break;
                               }
                               else continue;
                           }
                       }

                   }
                else continue;
            }
            else continue;
        }

        gsl_matrix *VV = gsl_matrix_calloc(p1a, p1a);
        gsl_matrix_memcpy(VV, sig);

        for(i=0;i<p1a;i++)
        {
            for(j=0;j<i;j++)    gsl_matrix_set(VV,i,j,gsl_matrix_get(VV,j,i));
        }
        gsl_permutation * vp = gsl_permutation_calloc(p1a);
        gsl_linalg_LU_decomp (VV, vp, &i);
        double u=sqrt(gsl_linalg_LU_det(VV,i));

        //Rprintf("det=%f\n",u);

        m=0;
        for(j=0;j<k;j++)
        {
            q=(int)gsl_vector_get(M1,j);

            cuh01=gsl_vector_get(CUH01,j);
            haz01=gsl_vector_get(HAZ01,j);

            for(i=0;i<p2;i++)
            {
                gsl_vector_set(X,i,gsl_matrix_get(C,j,2+i));
                gsl_vector_set(gammai,i,gsl_matrix_get(gamma,0,i));
            }
            xgamma1=MulVV(X,gammai);

            /*Extract Bayes estimate*/
            for (i=0;i<p1a;i++) gsl_vector_set(bi, i, gsl_matrix_get(Posbi, j, i));
            for (i=0;i<p1a;i++)
            {
                gsl_matrix_get_row(ci, Poscov, j*p1a+i);
                gsl_matrix_set_row(covi, i, ci);
            }
            temp1=0;

            if(p1a==1)
            {

                for(db0=0;db0<point;db0++)
                {
                        gsl_vector_set(rii,0,gsl_vector_get(xi, db0));
                        MulM(covi, rii, tii);
                        gsl_vector_scale(tii, sqrt(2));
                        gsl_vector_add(tii, bi);
                        gsl_vector_memcpy(ti, tii);

                        temp=1;

                        for(i=0;i<q;i++)
                        {
                            for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                            for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                            mu=MulVV(Z,beta);

                            temp*=1/sqrt(sigma*2*M_PI)*exp(-1/(2*sigma)*gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-MulVV(ti,xtilde)));
                        }

                        if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+MulVV(vee1,ti));

                        temp*=exp(0-cuh01*exp(xgamma1+MulVV(vee1,ti)));
                        temp*=1/sqrt(M_PI)*gsl_vector_get(wi,db0);
                        for (i=0;i<p1a;i++) temp*=gsl_matrix_get(covi, i, i);
                        temp/=u;
                        gsl_blas_dgemv(CblasNoTrans, 1.0, VC, tii, 0.0, tiii);
                        gsl_blas_ddot(tii, tiii, &uu);
                        temp*=exp(gsl_pow_2(gsl_blas_dnrm2(rii)) - uu/2);

                        temp1+=temp;


                }
            }


            if(p1a==2)
            {

                for(db0=0;db0<point;db0++)
                {
                    for(db1=0;db1<point;db1++)
                    {
                        gsl_vector_set(rii,0,gsl_vector_get(xi, db0));
                        gsl_vector_set(rii,1,gsl_vector_get(xi, db1));
                        MulM(covi, rii, tii);
                        gsl_vector_scale(tii, sqrt(2));
                        gsl_vector_add(tii, bi);
                        gsl_vector_memcpy(ti, tii);

                        temp=1;

                        for(i=0;i<q;i++)
                        {
                            for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                            for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                            mu=MulVV(Z,beta);

                        temp*=1/sqrt(sigma*2*M_PI)*exp(-1/(2*sigma)*gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-MulVV(ti,xtilde)));
                        }

                        if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+MulVV(vee1,ti));

                        temp*=exp(0-cuh01*exp(xgamma1+MulVV(vee1,ti)));
                        temp*=1/sqrt(gsl_pow_2(M_PI))*gsl_vector_get(wi,db0)*gsl_vector_get(wi,db1);
                        for (i=0;i<p1a;i++) temp*=gsl_matrix_get(covi, i, i);
                        temp/=u;
                        gsl_blas_dgemv(CblasNoTrans, 1.0, VC, tii, 0.0, tiii);
                        gsl_blas_ddot(tii, tiii, &uu);
                        temp*=exp(gsl_pow_2(gsl_blas_dnrm2(rii)) - uu/2);

                        temp1+=temp;

                    }
                }
            }

            if(p1a==3)
            {
                for(db0=0;db0<point;db0++)
                {
                    for(db1=0;db1<point;db1++)
                    {
                        for(db2=0;db2<point;db2++)
                        {

                            gsl_vector_set(rii,0,gsl_vector_get(xi, db0));
                            gsl_vector_set(rii,1,gsl_vector_get(xi, db1));
                            gsl_vector_set(rii,2,gsl_vector_get(xi, db2));
                            MulM(covi, rii, tii);
                            gsl_vector_scale(tii, sqrt(2));
                            gsl_vector_add(tii, bi);
                            gsl_vector_memcpy(ti, tii);

                            temp=1;

                            for(i=0;i<q;i++)
                            {
                                for(t=0;t<p1;t++)  gsl_vector_set(Z,t,gsl_matrix_get(Y,m+i,t+p1a+1));
                                for(t=0;t<p1a;t++)  gsl_vector_set(xtilde,t,gsl_matrix_get(Y,m+i,t+1));

                                mu=MulVV(Z,beta);

                            temp*=1/sqrt(sigma*2*M_PI)*exp(-1/(2*sigma)*gsl_pow_2(gsl_matrix_get(Y,m+i,0)-mu-MulVV(ti,xtilde)));
                            }

                            if(gsl_matrix_get(C,j,1)==1)  temp*=haz01*exp(xgamma1+MulVV(vee1,ti));

                            temp*=exp(0-cuh01*exp(xgamma1+MulVV(vee1,ti)));
                            temp*=1/sqrt(gsl_pow_3(M_PI))*gsl_vector_get(wi,db0)*gsl_vector_get(wi,db1)*gsl_vector_get(wi,db2);
                            for (i=0;i<p1a;i++) temp*=gsl_matrix_get(covi, i, i);
                            temp/=u;
                            gsl_blas_dgemv(CblasNoTrans, 1.0, VC, tii, 0.0, tiii);
                            gsl_blas_ddot(tii, tiii, &uu);
                            temp*=exp(gsl_pow_2(gsl_blas_dnrm2(rii)) - uu/2);

                            temp1+=temp;


                        }
                    }
                }
            }


            loglik+=log(temp1);

            m+=q;

        }


        gsl_vector_free(Z);
        gsl_vector_free(X);
        gsl_vector_free(xtilde);
        gsl_vector_free(gammai);

        gsl_vector_free(xi);
        gsl_vector_free(ti);
        gsl_vector_free(wi);
        gsl_vector_free(bi);
        gsl_vector_free(ci);
        gsl_matrix_free(covi);
        gsl_matrix_free(VC);
        gsl_matrix_free(VV);
        gsl_vector_free(rii);
        gsl_vector_free(tii);
        gsl_vector_free(tiii);
        gsl_vector_free(CUH01);
        gsl_vector_free(HAZ01);
        gsl_vector_free(CumuH01);
        gsl_permutation_free(vp);


        return loglik;
    }

    int Diff(
             const gsl_vector *prebeta,
             const gsl_vector *beta,
             const gsl_matrix *pregamma,
             const gsl_matrix *gamma,
             const gsl_vector *prevee1,
             const gsl_vector *vee1,
             const gsl_matrix *preH01,
             const gsl_matrix *H01,
             const double presigma,
             const double sigma,
             const gsl_matrix *presig,
             const gsl_matrix *sig,
             const double tol
             )
    {

         double epsilon=tol;

         if(DiffV(prebeta,beta)>epsilon || DiffM(pregamma,gamma)>epsilon || DiffV(prevee1,vee1)>epsilon
            || DiffM1(preH01,H01)==1 || Abs(presigma,sigma)>epsilon || DiffM(presig,sig)>epsilon)

         return 1;

         else return 0;

    }

    double DiffV(const gsl_vector *veca, const gsl_vector *vecb)
    {
        int k=veca->size;
        int i;
        double diff=0;

        for(i=0;i<k;i++)
        {
            if(Abs(gsl_vector_get(veca,i),gsl_vector_get(vecb,i))>diff) diff=Abs(gsl_vector_get(veca,i),gsl_vector_get(vecb,i));
        }

        return (diff);
    }

    double DiffM(const gsl_matrix *matrixa, const gsl_matrix *matrixb)
    {
        int nrow=matrixa->size1, ncol=matrixa->size2;
        int i, j;
        double diff=0;

        for(i=0;i<nrow;i++)
        {
            for(j=0;j<ncol;j++)
            {
                if(Abs(gsl_matrix_get(matrixa,i,j),gsl_matrix_get(matrixb,i,j))>diff)
                   diff=Abs(gsl_matrix_get(matrixa,i,j),gsl_matrix_get(matrixb,i,j));
            }
        }

        return (diff);
    }

    double Abs(const double a, const double b)
    {

        if (a>=b) return a-b;
        else return b-a;
    }

    void MulV(const gsl_vector *Z,gsl_matrix *ZZ)
    {
        int p = Z->size;
        int i,j;

        for(i=0;i<p;i++)
        {
            for(j=0;j<p;j++) gsl_matrix_set(ZZ,i,j,gsl_vector_get(Z,i)*gsl_vector_get(Z,j));
        }
    }


    double MulVV(const gsl_vector *Z,const gsl_vector *beta)
    {
        int p=Z->size;
        int i;
        double temp=0;

        for(i=0;i<p;i++)  temp+=gsl_vector_get(Z,i)*gsl_vector_get(beta,i);

        return (temp);
    }


    void MulMM(const gsl_matrix *A,const gsl_matrix *B,gsl_matrix *AB)
    {
        int p=A->size1;
        int q=A->size2;
        int k=B->size2;

        int i,j,t;
        double temp;

        for(i=0;i<p;i++)
        {
            for(j=0;j<k;j++)
            {
                temp=0;
                for(t=0;t<q;t++)  temp+=gsl_matrix_get(A,i,t)*gsl_matrix_get(B,t,j);
                gsl_matrix_set(AB,i,j,temp);
            }
        }

    }

    double Min(const double t1, const double t2)
    {
        if(t1<t2) return t1;
        else return t2;
    }

    int DiffM1(const gsl_matrix *matrixa, const gsl_matrix *matrixb)
    {
        int nrow=matrixa->size1, ncol=matrixa->size2;
        int i, j;
        int diff=0;

        i=nrow-1;

            for(j=0;j<ncol;j++)
            {
                if(gsl_matrix_get(matrixa,i,j)-gsl_matrix_get(matrixb,i,j)>0.001 || gsl_matrix_get(matrixa,i,j)-gsl_matrix_get(matrixb,i,j)<-0.001)

                    diff=1;

            }

        return diff;
    }

    void TransM(const gsl_matrix *A, gsl_matrix *B)
    {
        int rowa = A->size1;
        int cola = A->size2;

        int i, j;

        for(i=0;i<rowa;i++)
        {
            for(j=0;j<cola;j++)
            {
                gsl_matrix_set(B,j,i,gsl_matrix_get(A,i,j));
            }
        }

    }

    void STAT(gsl_matrix *store,int i,double *mean,double *sd)
    {
        int n=store->size1;
        int j;

        *mean=0, *sd=0;

        for(j=0;j<n;j++)   *mean+=gsl_matrix_get(store,j,i);
        *mean=*mean/(double)n;

        for(j=0;j<n;j++)   *sd+=gsl_pow_2(gsl_matrix_get(store,j,i)-*mean);
        *sd = sqrt(*sd/(double)(n-1));

    }

    int GetN(double t)
    {

        return (int)(t/0.5+1);

    }



    int inv_matrix(gsl_matrix *x_square)
    {
       int i,j;
       int k = x_square->size1;

       int status;

       gsl_vector *temp_vector=gsl_vector_calloc(k),
                  *solution=gsl_vector_calloc(k);
       gsl_matrix *out = gsl_matrix_calloc(k,k);

       for(i=0;i<k;i++)
       {
           for(j=0;j<k;j++) gsl_matrix_set(out,i,j,gsl_matrix_get(x_square,i,j));
       }

       status=gsl_linalg_cholesky_decompn(out);
       if(status==100) return status;

       for (i = 0; i < k; i++)
       {
           gsl_vector_set_all(temp_vector,0);
           gsl_vector_set(temp_vector,i,1);

           status=gsl_linalg_cholesky_solven(out, temp_vector, solution);
           if(status==100) return status;

           gsl_matrix_set_col(x_square,i,solution);
       }

       gsl_vector_free(temp_vector);
       gsl_vector_free(solution);
       gsl_matrix_free(out);

       return 0;

    }

    Rcpp::List jmcsf_cmain(double tol, int k, int n1,int p1,int p2, int p1a, int maxiter, int point,std::vector<double> xs,  std::vector<double> ws, std::string yfile, std::string cfile, std::string mfile, std::string Betasigmafile, std::string Sigcovfile,
                           std::vector<double> gammafile, int trace)
    {
        int g=1;
        /* allocate space for data */
        gsl_matrix *C = gsl_matrix_calloc(k,p2+2);
        gsl_matrix *Y= gsl_matrix_calloc(n1, p1+p1a+1);
        gsl_vector *M1= gsl_vector_calloc(k);
        gsl_matrix *Posbi=gsl_matrix_calloc(k, p1a);
        gsl_matrix *Poscov=gsl_matrix_calloc(k*p1a, p1a);
        gsl_matrix *Sigcov=gsl_matrix_calloc(p1a, p1a);
        gsl_vector *Betasigma=gsl_vector_calloc(p1+1);
        gsl_matrix * Cov=gsl_matrix_calloc(p1+g*p2+g*p1a+p1a*(p1a+1)/2+1,p1+g*p2+g*p1a+p1a*(p1a+1)/2+1);
        Rprintf("The survival dataset has single failure!\n");

        /* read Y matrix  */
         {
           FILE * f = fopen(yfile.c_str(), "r");

           if (f == NULL)
           {
               Rprintf("File %s does not exist.\n", yfile.c_str());
             return R_NilValue;
           }


           int nrows=0;
             // Extract characters from file and store in character c
           for (char c = fgetc(f); c != EOF; c = fgetc(f))
                   if (c == '\n')  nrows = nrows + 1;
           nrows=nrows+1;
           if (n1==nrows)
           {   rewind(f);
               gsl_matrix_fscanf(f, Y);
               fclose(f);
           }
           else
           {
               Rprintf("Input oberservations is %d, but the number of rows in %s is %d",n1,yfile.c_str(),nrows);
               fclose(f);
               return R_NilValue;
           }



         }

         /* read M1 vector  */
         {
             FILE * f = fopen(mfile.c_str(), "r");

             if (f == NULL)
             {
                 Rprintf("File %s does not exist.\n", mfile.c_str());
                 return R_NilValue;
             }


             int nrows=0;
             // Extract characters from file and store in character c
             for (char c = fgetc(f); c != EOF; c = fgetc(f))
                 if (c == '\n')  nrows = nrows + 1;
                 nrows=nrows+1;
                 if (k==nrows)
                 {   rewind(f);
                     gsl_vector_fscanf(f, M1);
                     fclose(f);
                 }
                 else
                 {
                     Rprintf("Input subjects is %d, but the number of rows in %s is %d",k,mfile.c_str(),nrows);
                     fclose(f);
                     return R_NilValue;
                 }


         }

         /* read Betasigma vector  */
         {
             FILE * f = fopen(Betasigmafile.c_str(), "r");

             if (f == NULL)
             {
                 Rprintf("File %s does not exist.\n", Betasigmafile.c_str());
                 return R_NilValue;
             }


             int nrows=0;
             // Extract characters from file and store in character c
             for (char c = fgetc(f); c != EOF; c = fgetc(f))
                 if (c == '\n')  nrows = nrows + 1;
                 nrows=nrows+1;
                 if (p1+1==nrows)
                 {   rewind(f);
                     gsl_vector_fscanf(f, Betasigma);
                     fclose(f);
                 }
                 else
                 {
                     Rprintf("Input subjects is %d, but the number of rows in %s is %d",p1+1,Betasigmafile.c_str(),nrows);
                     fclose(f);
                     return R_NilValue;
                 }


         }

         /* read Sigcov matrix  */
         {
             FILE * f = fopen(Sigcovfile.c_str(), "r");

             if (f == NULL)
             {
                 Rprintf("File %s does not exist.\n", Sigcovfile.c_str());
                 return R_NilValue;
             }


             int nrows=0;
             // Extract characters from file and store in character c
             for (char c = fgetc(f); c != EOF; c = fgetc(f))
                 if (c == '\n')  nrows = nrows + 1;
                 nrows=nrows+1;
                 if (p1a==nrows)
                 {   rewind(f);
                     gsl_matrix_fscanf(f, Sigcov);
                     fclose(f);
                 }
                 else
                 {
                     Rprintf("Input subjects is %d, but the number of rows in %s is %d",p1a,Sigcovfile.c_str(),nrows);
                     fclose(f);
                     return R_NilValue;
                 }


         }
         /* read C matrix  */
         {
           FILE * f = fopen(cfile.c_str(), "r");

           if (f == NULL)
           {
               Rprintf("File %s does not exist.\n", cfile.c_str());
             return R_NilValue;
           }

           int nrows=0;
             // Extract characters from file and store in character c
           for (char c = fgetc(f); c != EOF; c = fgetc(f))
                   if (c == '\n')  nrows = nrows + 1;
                   nrows=nrows+1;
                   if (k==nrows)
                    {   rewind(f);
                        gsl_matrix_fscanf(f, C);
                        fclose(f);
                    }
                    else
                    {
                        Rprintf("Input subjects is %d, but the number of rows in %s is %d",k,cfile.c_str(),nrows);
                        fclose(f);
                        return R_NilValue;
                    }


         }

        int i,j,iter,status;
        /* allocate space for estimated parameters */
        gsl_matrix * gamma=gsl_matrix_calloc(g, p2);
        gsl_vector * vee1=gsl_vector_calloc(p1a);
        gsl_vector * beta=gsl_vector_calloc(p1);
        gsl_matrix * sig=gsl_matrix_calloc(p1a, p1a);

        double sigma;
        for (i=0;i<p1;i++) gsl_vector_set(beta, i, gsl_vector_get(Betasigma, i));
        sigma = gsl_vector_get(Betasigma, p1);
        /*calculate the Posbi and Poscov*/
        i=0;
        //auto start_Pos = std::chrono::high_resolution_clock::now();
        gsl_matrix *Xinv = gsl_matrix_calloc(p1,p1);
        GetPosbi(Y,beta,M1,Sigcov,sigma,Xinv,Posbi,k,p1,p1a);
        inv_matrix(Xinv);
        GetPoscov(Y,beta,M1,Sigcov,sigma,Xinv,Poscov,k,p1,p1a);
        /* allocate space for pre parameters */

        gsl_vector * prebeta=gsl_vector_calloc(p1);
        gsl_matrix * pregamma=gsl_matrix_calloc(g,p2);
        gsl_vector * prevee1=gsl_vector_calloc(p1a);
        gsl_matrix * presig=gsl_matrix_calloc(p1a, p1a);

        double presigma;

        /*allocate space for standard error estimate*/
        gsl_vector * vbeta=gsl_vector_calloc(p1);
        gsl_matrix * vgamma=gsl_matrix_calloc(g,p2);
        gsl_vector * vvee1=gsl_vector_calloc(p1a);
        gsl_vector * vsig=gsl_vector_calloc(p1a*(p1a+1)/2);
        double v_sigma;

        //Define C_new matrix for risk set calculation
        gsl_matrix *C_new = gsl_matrix_calloc(k,2+p2);
        gsl_vector *p = gsl_vector_calloc(k);
        gsl_vector *pi = gsl_vector_calloc(2+p2);
        gsl_permutation * perm = gsl_permutation_calloc(k);
        gsl_permutation * rank = gsl_permutation_calloc(k);
        //return rank index
        //auto start_rank = std::chrono::high_resolution_clock::now();
        gsl_matrix_get_col(p, C, 0);
        gsl_sort_vector_index(perm, p);
        gsl_permutation_inverse(rank, perm);

        int index;
        for (j=0;j<k;j++)
        {
          gsl_matrix_get_row(pi, C, j);
          index = k - 1 - (int) gsl_permutation_get(rank, j);
          gsl_matrix_set_row(C_new, index, pi);
        }


        gsl_matrix * FH01 = gsl_matrix_calloc(2,k);                 /*** n = 800 here can be more flexible *****/

        gsl_matrix_set_zero(FH01);

        /* find # events for risk 1 */
        int u,a=0;
        u=0;
        for (j=0;j<k;j++)
        {
            if (gsl_matrix_get(C_new,j,1) == 1)
            {
                u++;
                if (j == k-1)
                {
                    a++;
                    gsl_matrix_set(FH01,0,k-a, gsl_matrix_get(C_new,j,0));
                    gsl_matrix_set(FH01,1,k-a,u);
                    u=0;
                }
                else if (gsl_matrix_get(C_new,j+1,0) != gsl_matrix_get(C_new,j,0))
                {
                    a++;
                    gsl_matrix_set(FH01,0,k-a, gsl_matrix_get(C_new,j,0));
                    gsl_matrix_set(FH01,1,k-a,u);
                    u=0;
                }
                else
                {
                    for (j=j+1;j<k;j++)
                    {
                        if (gsl_matrix_get(C_new,j,1) == 1)
                        {
                            u++;
                            if (j == k-1)
                            {
                                a++;
                                gsl_matrix_set(FH01,0,k-a, gsl_matrix_get(C_new,j,0));
                                gsl_matrix_set(FH01,1,k-a,u);
                                u=0;
                                break;
                            }
                            else if (gsl_matrix_get(C_new,j+1,0) != gsl_matrix_get(C_new,j,0))
                            {
                                a++;
                                gsl_matrix_set(FH01,0,k-a, gsl_matrix_get(C_new,j,0));
                                gsl_matrix_set(FH01,1,k-a,u);
                                u=0;
                                break;
                            }
                            else continue;
                        }
                        else
                        {
                            if (j == k-1)
                            {
                                a++;
                                gsl_matrix_set(FH01,0,k-a, gsl_matrix_get(C_new,j,0));
                                gsl_matrix_set(FH01,1,k-a,u);
                                u=0;
                                break;
                            }
                            else if (gsl_matrix_get(C_new,j+1,0) != gsl_matrix_get(C_new,j,0))
                            {
                                a++;
                                gsl_matrix_set(FH01,0,k-a, gsl_matrix_get(C_new,j,0));
                                gsl_matrix_set(FH01,1,k-a,u);
                                u=0;
                                break;
                            }
                            else continue;
                        }
                    }
                }

            }
            else continue;
        }

        if(a==0)
        {
            printf("No failure time information for risk 1; Program exits\n");
            return R_NilValue;
        }
        gsl_matrix * H01 = gsl_matrix_calloc(3,a);
        for(i=0;i<3;i++)
        {
            if(i<=1)
            {
                for(j=a;j>0;j--)    gsl_matrix_set(H01,i,a-j, gsl_matrix_get(FH01,i,k-j));
            }
            if(i==2)
            {
                for(j=0;j<a;j++)    gsl_matrix_set(H01,i,j,0.0001);
            }
        }

        gsl_matrix * preH01 = gsl_matrix_calloc(3,a);

        /* initialize the parameters */

        for(i=0;i<p2;i++)   gsl_matrix_set(gamma, 0, i, gammafile[i]);
        gsl_vector_set_zero(vee1);
        gsl_matrix_memcpy(sig, Sigcov);

        int index1;
        /*Sort data*/
        gsl_vector *M1_new = gsl_vector_calloc(k);
        gsl_vector *Posbivec = gsl_vector_calloc(p1a);
        gsl_matrix *Posbi_new = gsl_matrix_calloc(k,p1a);
        gsl_matrix *Poscov_new = gsl_matrix_calloc(k*p1a,p1a);
        for (j=0;j<k;j++)
        {
          gsl_matrix_get_row(Posbivec,Posbi,j);
          index1 = k - 1 - (int) gsl_permutation_get(rank, j);
          gsl_vector_set(M1_new,index1,gsl_vector_get(M1,j));
          gsl_matrix_set_row(Posbi_new, index1, Posbivec);
        }
        //Define rank of long data for each subject
        gsl_vector *M1_rank = gsl_vector_calloc(k);
        gsl_vector *M1_new_rank = gsl_vector_calloc(k);
        int n2 = Y->size2;
        gsl_vector *Y_row = gsl_vector_calloc(n2);
        gsl_matrix *Y_new = gsl_matrix_calloc(n1,n2);
        double temp1=0;
        for (j=1;j<k;j++)
        {
            temp1+=gsl_vector_get(M1,j-1);
            gsl_vector_set(M1_rank,j,temp1);
        }
        temp1=0;
        for (j=1;j<k;j++)
        {
            temp1+=gsl_vector_get(M1_new,j-1);
            gsl_vector_set(M1_new_rank,j,temp1);
        }

        int u1,q;
        for (j=0;j<k;j++)
        {
            index1 = k - 1 - (int) gsl_permutation_get(rank, j);
            u1 = gsl_vector_get(M1,j);
            for (q=0;q<u1;q++)
            {
                gsl_matrix_get_row(Y_row, Y, gsl_vector_get(M1_rank, j)+q);
                gsl_matrix_set_row(Y_new, gsl_vector_get(M1_new_rank, index1)+q, Y_row);
            }
        }

        for (j=0;j<k;j++)
        {
            index1 = k - 1 - (int) gsl_permutation_get(rank, j);
            for (q=0;q<p1a;q++)
            {
                gsl_matrix_get_row(Posbivec,Poscov,j*p1a+q);
                gsl_matrix_set_row(Poscov_new, index1*p1a+q, Posbivec);
            }

        }

        gsl_matrix_free(Y);
        gsl_matrix_free(C);
        gsl_vector_free(M1);
        gsl_matrix_free(Poscov);
        gsl_matrix_free(Posbi);
        gsl_vector_free(Posbivec);
        gsl_vector_free(Y_row);

        iter=0;
        do
        {
            iter=iter+1;

            /* store the pre-information */

            gsl_vector_memcpy(prebeta, beta);
            gsl_vector_memcpy(prevee1, vee1);
            gsl_matrix_memcpy(pregamma, gamma);
            gsl_matrix_memcpy(preH01, H01);
            gsl_matrix_memcpy(presig, sig);
            presigma=sigma;

            if (trace == 1 && iter == 1)
            {
                Rprintf("iter=%d   status=%d\n",iter,status);
                Rprintf("Beta = \n");
                for (i=0;i<p1;i++)
                {
                    Rprintf("%f     ", gsl_vector_get(beta,i));
                }
                Rprintf("\n");

                Rprintf("Gamma = \n");
                for (i=0;i<g;i++)
                {
                    for(j=0;j<p2;j++)
                    {
                        Rprintf("%f    ", gsl_matrix_get(gamma,i,j));
                    }
                    Rprintf("\n");
                }

                Rprintf("Vee1 = \n");
                for (i=0;i<p1a;i++)
                {
                    Rprintf("%f     ", gsl_vector_get(vee1,i));
                }
                Rprintf("\n");

                Rprintf("sigma = %f\n",sigma);

                Rprintf("Sig = \n");
                for (i=0;i<p1a;i++)
                {
                    for(j=0;j<p1a;j++)
                    {
                        Rprintf("%f    ", gsl_matrix_get(sig,i,j));
                    }
                    Rprintf("\n");
                }
            }

            /* get new parameter estimates */
            //auto start_EMstep = std::chrono::high_resolution_clock::now();

            status = EM(beta,gamma,vee1,H01,&sigma,sig,p1a,Y_new,C_new,M1_new,Posbi_new,Poscov_new,point,xs,ws);

            if (trace == 1)
            {
                Rprintf("iter=%d   status=%d\n",iter,status);
                Rprintf("Beta = \n");
                for (i=0;i<p1;i++)
                {
                    Rprintf("%f     ", gsl_vector_get(beta,i));
                }
                Rprintf("\n");

                Rprintf("Gamma = \n");
                for (i=0;i<g;i++)
                {
                    for(j=0;j<p2;j++)
                    {
                        Rprintf("%f    ", gsl_matrix_get(gamma,i,j));
                    }
                    Rprintf("\n");
                }

                Rprintf("Vee1 = \n");
                for (i=0;i<p1a;i++)
                {
                    Rprintf("%f     ", gsl_vector_get(vee1,i));
                }
                Rprintf("\n");

                Rprintf("sigma = %f\n",sigma);

                Rprintf("Sig = \n");
                for (i=0;i<p1a;i++)
                {
                    for(j=0;j<p1a;j++)
                    {
                        Rprintf("%f    ", gsl_matrix_get(sig,i,j));
                    }
                    Rprintf("\n");
                }
            }

        }while(Diff(prebeta,beta,pregamma,gamma,prevee1,vee1,preH01,H01,presigma,sigma,
               presig,sig, tol)==1
               && status != 100 && iter<maxiter);

        if(status==100)
        {
            Rprintf("program stops because of error\n");
            return R_NilValue;
        }
        if(iter==maxiter)
        {
            printf("program stops because of nonconvergence\n");
            return R_NilValue;
        }

        NumericMatrix vcmatrix(Cov->size1,Cov->size1);
        NumericVector betas(p1);
        NumericVector se_betas(p1);
        NumericMatrix gamma_matrix(g,p2);
        NumericMatrix sd_gamma_matrix(g,p2);
        NumericVector vee1_estimate(p1a);
        NumericVector sd_vee1_estimate(p1a);
        double sigma2_val;
        double se_sigma2_val;
        int iter_val;
        NumericMatrix sigma_matrix(p1a,p1a);
        NumericVector sd_sigma((p1a)*(p1a+1)/2);

        double loglike=0.0;
        if(status != 100 && iter<maxiter)
        {
            /* if algorithm coverges, compute the variance-covariance matrix of parameters ***/
            status = GetCov(Cov,beta,gamma,vee1,H01,sigma,sig,Y_new,C_new,M1_new,Posbi_new,Poscov_new,p1a,point,xs,ws);

            if(status==100)
            {
                printf("program stops because of error\n");
                return R_NilValue;
            }

            if(status != 100)
            {

                status=inv_matrix(Cov);

                if(status==100)
                {
                    printf("program stops because of error\n");
                    return R_NilValue;
                }

                if(status!=100)
                {
                    for (i=0;(unsigned)i<Cov->size1;i++)
                    {
                        for(j=0;(unsigned)j<Cov->size1;j++)
                        {
                            vcmatrix(i,j)=gsl_matrix_get(Cov,i,j);
                        }
                    }

                for (i=0;i<p1;i++)  gsl_vector_set(vbeta,i,gsl_matrix_get(Cov,i,i));
                v_sigma=gsl_matrix_get(Cov,p1,p1);

                for (i=p1+1;i<p1+p2+1;i++)  gsl_matrix_set(vgamma,0,i-p1-1,gsl_matrix_get(Cov,i,i));
                for (i=p1+p2+1;i<p1+p2+p1a+1;i++) gsl_vector_set(vvee1, i-p1-p2-1, gsl_matrix_get(Cov,i,i));
                for (i=p1+p2+p1a+1;i<p1+p2+p1a+1+p1a*(p1a+1)/2;i++) gsl_vector_set(vsig, i-p1-p2-p1a-1, gsl_matrix_get(Cov, i,i));

                for (i=0;i<p1;i++)
                {
                    betas(i)=gsl_vector_get(beta,i);
                }
                for (i=0;i<p1;i++)
                {
                    se_betas(i)=sqrt(gsl_vector_get(vbeta,i));
                }

                loglike = Getloglik(beta,gamma,vee1,H01,sigma,sig,Y_new,C_new,M1_new,Posbi_new,Poscov_new,p1a,point,xs,ws);

                for (i=0;i<g;i++)
                {
                    for(j=0;j<p2;j++)
                    {
                          gamma_matrix(i,j)=gsl_matrix_get(gamma,i,j);
                    }
                }

                for (i=0;i<g;i++)
                {
                   for(j=0;j<p2;j++)
                   {
                       sd_gamma_matrix(i,j)=sqrt(gsl_matrix_get(vgamma,i,j));
                   }
                }

                for (i=0;i<p1a;i++)
                {
                    vee1_estimate(i)=gsl_vector_get(vee1,i);
                }
                for (i=0;i<p1a;i++)
                {
                    sd_vee1_estimate(i)=sqrt(gsl_vector_get(vvee1,i));
                }

                sigma2_val=sigma;
                se_sigma2_val=sqrt(v_sigma);

                for (i=0;i<p1a;i++)
                {
                    for (j=0;j<p1a;j++)
                    {
                        sigma_matrix(i,j)=gsl_matrix_get(sig,i,j);
                    }
                }

                for (i=0;i<(p1a)*(p1a+1)/2;i++)
                {
                    sd_sigma(i)=sqrt(gsl_vector_get(vsig,i));
                }
                iter_val = iter;

                }
            }
        }

        gsl_matrix_free(FH01);
        gsl_matrix_free(H01);
        gsl_matrix_free(preH01);

        gsl_matrix_free(Y_new);
        gsl_vector_free(M1_new);

        gsl_matrix_free(gamma);
        gsl_vector_free(beta);
        gsl_vector_free(vee1);

        gsl_matrix_free(pregamma);
        gsl_vector_free(prebeta);
        gsl_vector_free(prevee1);

        gsl_matrix_free(vgamma);
        gsl_vector_free(vbeta);
        gsl_vector_free(vvee1);
        gsl_matrix_free(sig);
        gsl_matrix_free(presig);
        gsl_vector_free(vsig);
        gsl_matrix_free(Cov);
        //free C_new matrix
        gsl_matrix_free(C_new);
        gsl_vector_free(p);
        gsl_vector_free(pi);
        gsl_permutation_free(perm);
        gsl_permutation_free(rank);

        Rcpp::List ret;
        ret["vcmatrix"] = vcmatrix;
        ret["betas"] = betas;
        ret["se_betas"] = se_betas;
        ret["gamma_matrix"] = gamma_matrix;
        ret["se_gamma_matrix"] = sd_gamma_matrix;
        ret["vee1_estimate"] = vee1_estimate;
        ret["sd_vee1_estimate"] = sd_vee1_estimate;
        ret["sigma2_val"] = sigma2_val;
        ret["se_sigma2_val"] = se_sigma2_val;
        ret["sigma_matrix"] = sigma_matrix;
        ret["se_sigma"] = sd_sigma;
        ret["loglike"] = loglike;
        ret["iters"] = iter_val;
        return ret;

    }




}
