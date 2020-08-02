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

void MulM(const gsl_matrix *XX, const gsl_vector *X, gsl_vector *beta);

int GetN(double t);

double Min(const double t1, const double t2);

int gsl_linalg_cholesky_svx (const gsl_matrix * LLT,
                             gsl_vector * x);

int gsl_linalg_cholesky_solven (const gsl_matrix * LLT,
                                const gsl_vector * b,
                                gsl_vector * x);

int gsl_linalg_cholesky_decompn (gsl_matrix * A);

