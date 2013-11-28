#ifndef __LA_H
#define __LA_H

#include <gsl/gsl_blas.h>

#define gsl_matrix_mult(A,B,C) (gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, B, 0, C))
#define LOGDET1M_TAYLOR_EPS 1e-12

double la_norm_froebenius(gsl_matrix *M);
double la_matrix_trace(gsl_matrix *A);


#endif
