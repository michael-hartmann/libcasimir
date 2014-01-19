#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "libcasimir.h"
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

/* get matrix element m_(i,j) */
__float inline matrix_get(const matrix_t *m, size_t i, size_t j)
{
    return m->M[i*m->size+j];
}

/* set matrix element m_(i,j) to x */
void inline matrix_set(matrix_t *m, size_t i, size_t j, __float x)
{
    m->M[i*m->size+j] = x;
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
                    s = 0.0q;
                }
                else if(a == 0)
                {
                    c = 0.0q;
                    s = -__copysign(1, b);
                }
                else if(__abs(b) > __abs(a))
                {
                    __float t = a/b;
                    __float u = __copysign(__sqrt(1+t*t),b);
                    s = -1.0q/u;
                    c = -s*t;
                }
                else
                {
                    __float t = b/a;
                    __float u = __copysign(__sqrt(1+t*t),a);
                    c = 1.0q/u;
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
