#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>

#include "la.h"

double la_norm_froebenius(gsl_matrix *M)
{
    int dim = M->size1;
    double norm = 0;
    int i,j;

    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
            norm += gsl_pow_2(gsl_matrix_get(M, i, j));

    return sqrt(norm);
}

double la_matrix_trace(gsl_matrix *A)
{
    double trace = 0;
    int i;

    for(i = 0; i < A->size1; i++)
        trace += gsl_matrix_get(A, i, i);

    return trace;
}

void la_matrix_info(FILE *stream, gsl_matrix *M)
{
    int i,j;
    int zeros = 0, nans = 0, infs = 0;
    double min, max, diag_min, diag_max;
    min = max = fabs(gsl_matrix_get(M, 0, 0));
    diag_min = min;
    diag_max = max;

    for(i = 0; i < M->size1; i++)
    {
        double elem = fabs(gsl_matrix_get(M, i,i));
        diag_min = MIN(diag_min, elem);
        diag_max = MAX(diag_max, elem);
    }

    for(i = 0; i < M->size1; i++)
        for(j = 0; j < M->size2; j++)
         {
             double elem = fabs(gsl_matrix_get(M, i,j));
             if(isnan(elem))
                 nans++;
             if(isinf(elem))
                 infs++;
             if(elem == 0)
                 zeros++;
             min = MIN(min, elem);
             max = MAX(max, elem);
         }

    //printf("diag_diff=%g\n", log(fabs(diag_max))-log(fabs(diag_min)));
    if(zeros || nans || infs)
        fprintf(stream, "zeros=%d, nans=%d, infs=%d, min=%g, max=%g, diag_min=%g, diag_max=%g, diff=%g\n", zeros, nans, infs, min, max, diag_min, diag_max, log10(fabs(diag_max))-log10(fabs(diag_min)));
}

