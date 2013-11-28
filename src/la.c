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
