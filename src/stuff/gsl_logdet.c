#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <math.h>

#define FLOAT_RADIX       2.0
#define FLOAT_RADIX_SQ    (FLOAT_RADIX * FLOAT_RADIX)
#define LOG_FLOAT_RADIX   M_LN2
#define LOG_FLOAT_RADIX_SQ 2*M_LN2
#define LOG_095 -0.05129329438755058

#include "gsl_logdet.h"

void gsl_log_balance(gsl_matrix *A)
{
    size_t i,j;
    const size_t N = A->size1;
    int not_converged = 1;
    double *D;
    double *list_row;
    double *list_column;

    D           = (double *)malloc(N*sizeof(double));
    list_row    = (double *)malloc(N*sizeof(double));
    list_column = (double *)malloc(N*sizeof(double));

    /* initialize D to the identity matrix */
    for(i = 0; i < N; i++)
        D[i] = 0;

    while(not_converged)
    {
        size_t len = 0;
        double g, f, s;
        double row_norm, col_norm;

        not_converged = 0; \

        for (i = 0; i < N; ++i)
        {
            len = 0;

            for (j = 0; j < N; ++j)
                if (j != i)
                {
                    list_column[len] = gsl_matrix_get(A,j,i);
                    list_row[len]    = gsl_matrix_get(A,i,j);
                    len++;
                }

            row_norm = logadd_m(list_row,    len);
            col_norm = logadd_m(list_column, len);

            if ((col_norm == log(0)) || (row_norm == log(0)))
              continue;

            g = row_norm - LOG_FLOAT_RADIX;
            f = 0;
            s = logadd(col_norm, row_norm);

            /* \
             * find the integer power of the machine radix which
             * comes closest to balancing the matrix
             */
            while (col_norm < g)
            {
                f += LOG_FLOAT_RADIX;
                col_norm += LOG_FLOAT_RADIX_SQ;
            }

            g = row_norm + LOG_FLOAT_RADIX;

            while (col_norm > g)
            {
                f -= LOG_FLOAT_RADIX;
                col_norm -= LOG_FLOAT_RADIX_SQ;
            }

            if (logadd(row_norm, col_norm) < (LOG_095+s+f))
            {
                int k;
                not_converged = 1;

                g = -f;

                /* \
                 * apply similarity transformation D, where
                 * D_{ij} = f_i * delta_{ij}
                 */

                /* multiply by D^{-1} on the left */
                for(k = 0; k < N; k++)
                    gsl_matrix_set(A, i,k, g+gsl_matrix_get(A,i,k));


                /* multiply by D on the right */
                for(k = 0; k < N; k++)
                    gsl_matrix_set(A, k,i, f+gsl_matrix_get(A,k,i));

                /* keep track of transformation */
                for(k = 0; k < N; k++)
                    D[k] += f;
            }
        }
    }

    free(D);
    free(list_column);
    free(list_row);
}

double gsl_logdet(gsl_matrix *M)
{
    size_t i;
    double value = 0;
    const int dim = M->size1;
    gsl_vector *tau = gsl_vector_alloc(dim);

    gsl_linalg_QR_decomp(M, tau);

    for(i = 0; i < dim; i++)
        value += log(fabs(gsl_matrix_get(M, i,i)));

    gsl_vector_free(tau);

    return value;
}
