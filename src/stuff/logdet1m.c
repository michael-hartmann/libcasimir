#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <math.h>

#include "logdet1m.h"


/*
* Calculate logarithm of determinant of E-M
* log(det(E-M)) = log((1-λ1)*(1-λ2)*...*(1-λn)) = log(1-λ1)+log(1-λ2)+...+log(1-λn)
* = Sum_i logm(-λi)
*
* where E is the identity matrix and λi are the eigenvalues of M
*/
int logdet1m_eigenvalues(gsl_matrix *M, double *logdet)
{
    int i, ret, dim = M->size1;
    *logdet = 0;

    gsl_vector_complex *eval = gsl_vector_complex_alloc(dim);

    gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc(dim);
    gsl_eigen_nonsymm_params(0, 1, w);

    if((ret = gsl_eigen_nonsymm(M, eval, w)) == 0)
    {
        for(i = 0; i < dim; i++)
        {
            gsl_complex lambda_i = gsl_vector_complex_get(eval, i);
            double real = GSL_REAL(lambda_i);
            double imag = GSL_IMAG(lambda_i);

            *logdet += 0.5*gsl_log1p( -2*real+gsl_pow_2(real) + gsl_pow_2(imag) );
        }
    }

    gsl_eigen_nonsymm_free(w);
    gsl_vector_complex_free(eval);

    return ret;
}
