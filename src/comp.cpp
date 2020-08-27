//
//  SimData.cpp
//  OldJMmodel
//
//  Created by Shanpeng Li on 6/22/20.
//
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
#include "comp.h"

#include <string>
#include <iostream>
#include <fstream>


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

