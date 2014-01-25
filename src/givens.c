#include <stdio.h>
#include <stdlib.h>

#include "libcasimir.h"
#include "sfunc.h"
#include "givens.h"

/* Allocate space for a quadratic matrix of size x size. */
matrix_t *matrix_alloc(size_t size)
{
    matrix_t *matrix = malloc(sizeof(matrix_t));
    if(matrix == NULL)
        return NULL;

    matrix->size = size;
    matrix->M = malloc(size*size*sizeof(__float));
    if(matrix->M == NULL)
    {
        matrix_free(matrix);
        return NULL;
    }

    return matrix;
}

/* Allocate space for a quadratic matrix of size x size. */
matrix_sign_t *matrix_sign_alloc(size_t size)
{
    matrix_sign_t *matrix = malloc(sizeof(matrix_sign_t));
    if(matrix == NULL)
        return NULL;

    matrix->size = size;
    matrix->M = malloc(size*size*sizeof(char));
    if(matrix->M == NULL)
    {
        matrix_sign_free(matrix);
        return NULL;
    }

    return matrix;
}

/*
void matrix_fprintf(const matrix_t *m, FILE *stream, const char *format, const char *lines, const char *rows)
{
    size_t i, j;
    size_t dim = m->size;

    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            fprintf(stream, format, matrix_get(m, i,j));
            fprintf(stream, "%s", rows);
        }
        fprintf(stream, "%s", lines);
    }
}
*/

/* Free space of matrix m */
void matrix_free(matrix_t *m)
{
    if(m->M != NULL)
    {
        free(m->M);
        m->M = NULL;
    }
    free(m);
}

/* Free space of matrix m */
void matrix_sign_free(matrix_sign_t *m)
{
    if(m->M != NULL)
    {
        free(m->M);
        m->M = NULL;
    }
    free(m);
}

/* calculate froebenius norm of matrix */
__float matrix_froebenius(matrix_t *M)
{
    int i,j;
    int dim = M->size;
    __float norm = 0;

    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
            norm += pow_2(matrix_get(M, i,j));

    return __sqrt(norm);
}

/* calculate log(det(M)) */
__float matrix_logdet(matrix_t *M)
{
    size_t i, j, n;
    size_t dim = M->size;
    __float det = 0;

    for(j = 0; j < dim-1; j++)
        for(i = j+1; i < dim; i++)
        {
            __float c, s;
            __float Mij = matrix_get(M, i, j);

            if(Mij != 0)
            {
                __float a = matrix_get(M, j, j);
                __float b = Mij;

                if(b == 0)
                {
                    c = __copysign(1,a);
                    s = 0;
                }
                else if(a == 0)
                {
                    c = 0;
                    s = -__copysign(1, b);
                }
                else if(__abs(b) > __abs(a))
                {
                    __float t = a/b;
                    __float u = __copysign(__sqrt(1+t*t),b);
                    s = -1/u;
                    c = -s*t;
                }
                else
                {
                    __float t = b/a;
                    __float u = __copysign(__sqrt(1+t*t),a);
                    c = 1/u;
                    s = -c*t;
                }

                for(n = 0; n < dim; n++)
                {
                    __float Min = matrix_get(M, i, n);
                    __float Mjn = matrix_get(M, j, n);

                    matrix_set(M, i, n, c*Min + s*Mjn);
                    matrix_set(M, j, n, c*Mjn - s*Min);
                }

                matrix_set(M, i, j, 0);
            }
        }

    for(i = 0; i < dim; i++)
        det += __log(__abs(matrix_get(M, i, i)));
    return det;
}

/*
void matrix_info(FILE *stream, matrix_t *M)
{
    int i,j;
    int zeros = 0, nans = 0, infs = 0;
    __float min, max, diag_min, diag_max;
    min = max = __abs(matrix_get(M, 0, 0));
    diag_min = min;
    diag_max = max;

    for(i = 0; i < M->size; i++)
    {
        __float elem = __abs(matrix_get(M, i,i));
        diag_min = MIN(diag_min, elem);
        diag_max = MAX(diag_max, elem);
    }

    for(i = 0; i < M->size; i++)
        for(j = 0; j < M->size; j++)
         {
             __float elem = __abs(matrix_get(M, i,j));
             if(isnan(elem))
                 nans++;
             if(isinf(elem))
                 infs++;
             if(elem == 0)
                 zeros++;
             min = MIN(min, elem);
             max = MAX(max, elem);
         }

    char str_diag_min[128];
    char str_diag_max[128];
    char str_max[120];
    char str_min[120];
    char str_diff_diag[120];
    quadmath_snprintf(str_diag_min,  sizeof(str_diag_min),  "%+Qg", diag_min);
    quadmath_snprintf(str_diag_max,  sizeof(str_diag_max),  "%+Qg", diag_max);
    quadmath_snprintf(str_min,       sizeof(str_min),       "%+Qg", min);
    quadmath_snprintf(str_max,       sizeof(str_max),       "%+Qg", max);
    quadmath_snprintf(str_diff_diag, sizeof(str_diff_diag), "%+Qg", (__log(diag_max)-__log(diag_min))/__log(10));
    fprintf(stream, "zeros=%d, nans=%d, infs=%d, min=%s, max=%s, diag_min=%s, diag_max=%s, diff=%s\n", zeros, nans, infs, str_min, str_max, str_diag_min, str_diag_max, str_diff_diag);
}
*/

__float matrix_absmax(matrix_t *M)
{
    int i,j;
    int dim = M->size;
    __float max = __abs(matrix_get(M, 0,0));

    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
            max = MAX(max, __abs(matrix_get(M, i,j)));

    return max;
}

__float matrix_absmin(matrix_t *M)
{
    int i,j;
    int dim = M->size;
    __float min = __abs(matrix_get(M, 0,0));

    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
            min = MIN(min, __abs(matrix_get(M, i,j)));

    return min;
}

/* Balance a general matrix by scaling the rows and columns, so the
 * new row and column norms are the same order of magnitude.
 *
 * B =  D^-1 A D
 *
 * where D is a diagonal matrix
 * 
 * This is necessary for the unsymmetric eigenvalue problem since the
 * calculation can become numerically unstable for unbalanced
 * matrices.  
 *
 * See Golub & Van Loan, "Matrix Computations" (3rd ed), Section 7.5.7
 * and Wilkinson & Reinsch, "Handbook for Automatic Computation", II/11 p320.
 */
void matrix_balance(matrix_t *A)
{
    size_t i,j;
    const size_t N = A->size;
    __float D[N];
    int not_converged = 1;

    /* initialize D to the identity matrix */
    for(i = 0; i < N; i++)
        D[i] = 1;

    while(not_converged)
    {
        __float g, f, s;
        __float row_norm, col_norm;

        not_converged = 0;

        for (i = 0; i < N; ++i)
        {
            row_norm = 0;
            col_norm = 0;

            for (j = 0; j < N; ++j)
                if (j != i)
                {
                  col_norm += __abs(matrix_get(A, j, i));
                  row_norm += __abs(matrix_get(A, i, j));
                }

            if ((col_norm == 0.0) || (row_norm == 0.0))
              continue;

            g = row_norm / FLOAT_RADIX;
            f = 1.0;
            s = col_norm + row_norm;

            /*
             * find the integer power of the machine radix which
             * comes closest to balancing the matrix
             */
            while (col_norm < g)
            {
                f *= FLOAT_RADIX;
                col_norm *= FLOAT_RADIX_SQ;
            }

            g = row_norm * FLOAT_RADIX;

            while (col_norm > g)
            {
                f /= FLOAT_RADIX;
                col_norm /= FLOAT_RADIX_SQ;
            }

            if ((row_norm + col_norm) < 0.95 * s * f)
            {
                int k;
                not_converged = 1;

                g = 1.0 / f;

                /*
                 * apply similarity transformation D, where
                 * D_{ij} = f_i * delta_{ij}
                 */

                /* multiply by D^{-1} on the left */
                for(k = 0; k < N; k++)
                    matrix_set(A, i,k, g*matrix_get(A,i,k));


                /* multiply by D on the right */
                for(k = 0; k < N; k++)
                    matrix_set(A, k,i, f*matrix_get(A,k,i));

                /* keep track of transformation */
                for(k = 0; k < N; k++)
                    D[k] *= f;
            }
        }
    }
}

void matrix_log_balance(matrix_t *A)
{
    size_t i,j;
    const size_t N = A->size;
    __float D[N];
    int not_converged = 1;

    /* initialize D to the identity matrix */
    for(i = 0; i < N; i++)
        D[i] = 0;

    while(not_converged)
    {
        __float g, f, s;
        __float row_norm, col_norm;

        not_converged = 0;

        for (i = 0; i < N; ++i)
        {
            row_norm = __log(0);
            col_norm = __log(0);

            for (j = 0; j < N; ++j)
                if (j != i)
                {
                    col_norm = logadd(col_norm, matrix_get(A,j,i));
                    row_norm = logadd(row_norm, matrix_get(A,i,j));
                }

            if ((col_norm == __log(0)) || (row_norm == __log(0)))
              continue;

            g = row_norm - LOG_FLOAT_RADIX;
            f = 0;
            s = logadd(col_norm, row_norm);

            /*
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

                /*
                 * apply similarity transformation D, where
                 * D_{ij} = f_i * delta_{ij}
                 */

                /* multiply by D^{-1} on the left */
                for(k = 0; k < N; k++)
                    matrix_set(A, i,k, g+matrix_get(A,i,k));


                /* multiply by D on the right */
                for(k = 0; k < N; k++)
                    matrix_set(A, k,i, f+matrix_get(A,k,i));

                /* keep track of transformation */
                for(k = 0; k < N; k++)
                    D[k] += f;
            }
        }
    }
}

void matrix_exp(matrix_t *m)
{
    size_t i, dim = m->size;
    __float *M = m->M;

    for(i = 0; i < dim*dim; i++)
        M[i] = __exp(M[i]);
}
