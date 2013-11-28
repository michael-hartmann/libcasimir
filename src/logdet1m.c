#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <math.h>

#include "la.h"
#include "logdet1m.h"


/*
int logdet1m_qr(gsl_matrix *M, double *logdet)
{
    int ret, i;
    int dim = M->size1;
    gsl_vector *tau = gsl_vector_alloc(dim);

    gsl_matrix_scale(M, -1);
    for(i = 0; i < dim; i++)
        gsl_matrix_set(M, i, i, 1+gsl_matrix_get(M, i, i));

    gsl_linalg_balance_matrix(M, tau);
    ret = gsl_linalg_QR_decomp(M, tau);
    gsl_vector_free(tau);

    if(ret == 0)
    {
        *logdet = 0;
        for(i = 0; i < dim; i++)
        {
            double elem = fabs(gsl_matrix_get(M,i,i));
            if(fabs(elem-1) < 1e-3)
                *logdet += gsl_log1p(elem-1);
            else
                *logdet += log(elem);
        }
    }

    return ret;
}
*/

int logdet1m_taylor(gsl_matrix *M, double *logdet)
{
    double logdet_last;
    int j, dim = M->size1, ret = 1;

    gsl_matrix *M_temp = gsl_matrix_alloc(dim,dim);
    gsl_matrix *M_pow  = gsl_matrix_alloc(dim,dim);
    gsl_matrix_memcpy(M_pow, M);

    *logdet = logdet_last = -la_matrix_trace(M_pow);

    for(j = 2; j < 250; j++)
    {
        gsl_matrix_mult(M_pow,M,M_temp);
        gsl_matrix_memcpy(M_pow, M_temp);

        *logdet -= la_matrix_trace(M_pow)/j;

        if(*logdet == 0 && logdet_last == 0)
        {
            ret = 0;
            break;
        }
        if(fabs(*logdet/logdet_last-1) < LOGDET1M_TAYLOR_EPS)
        {
            ret = 0;
            break;
        }

        logdet_last = *logdet;
    }

    gsl_matrix_free(M_temp);
    gsl_matrix_free(M_pow);

    return ret;
}

/*
int logdet_qr(gsl_matrix *M, double *logdet)
{
    int ret, i, dim = M->size1;
    gsl_vector *tau = gsl_vector_alloc(dim);

    gsl_linalg_balance_matrix(M, tau);
    ret = gsl_linalg_QR_decomp(M, tau);
    gsl_vector_free(tau);

    if(ret == 0)
    {
        *logdet = 0;
        for(i = 0; i < dim; i++)
        {
            double elem = fabs(gsl_matrix_get(M,i,i));
            //printf("%d %g\n", i, elem);
            if(fabs(elem-1) < 1e-3)
                *logdet += gsl_log1p(elem-1);
            else
                *logdet += log(elem);
        }
    }

    return ret;
}
*/

/*
 * Calculate logarithm of determinant of E-M
 * log(det(E-M)) = log((1-λ1)*(1-λ2)*...*(1-λn)) = log(1-λ1)+log(1-λ2)+...+log(1-λn)
 *               = Sum_i logm(-λi)
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
