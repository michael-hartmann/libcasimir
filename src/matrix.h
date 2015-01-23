#ifndef __MATRIX_H
#define __MATRIX_H

#include <stdio.h>
#include <math.h>

#include "edouble.h"
#include "utils.h"

#define FLOAT_RADIX       2.0
#define FLOAT_RADIX_SQ    (FLOAT_RADIX * FLOAT_RADIX)
#define LOG_FLOAT_RADIX   M_LN2
#define LOG_FLOAT_RADIX_SQ 2*M_LN2
#define LOG_095 -0.05129329438755058


#define MATRIX_TYPEDEF(NAME, MATRIX_TYPE) \
    typedef struct { \
        size_t size; \
        MATRIX_TYPE *M; \
    } NAME


MATRIX_TYPEDEF(matrix_t, double);
MATRIX_TYPEDEF(matrix_char_t, char);
MATRIX_TYPEDEF(matrix_edouble_t, edouble);

#define MATRIX_ALLOC(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) \
    MATRIX_TYPE *FUNCTION_PREFIX ## _alloc(size_t size)  \
    { \
        MATRIX_TYPE *matrix = xmalloc(sizeof(MATRIX_TYPE)); \
        if(matrix == NULL) \
            return NULL; \
 \
        matrix->size = size; \
        matrix->M = xmalloc(size*size*sizeof(TYPE)); \
        if(matrix->M == NULL) \
        { \
            FUNCTION_PREFIX ## _free(matrix); \
            return NULL; \
        } \
 \
        return matrix; \
    }

#define MATRIX_ALLOC_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) MATRIX_TYPE *FUNCTION_PREFIX ## _alloc(size_t size)

#define MATRIX_FREE(FUNCTION_PREFIX, MATRIX_TYPE) \
    void FUNCTION_PREFIX ## _free(MATRIX_TYPE *m) \
    { \
        if(m->M != NULL) \
        { \
            xfree(m->M); \
            m->M = NULL; \
        } \
        xfree(m); \
    }

#define MATRIX_FREE_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) void FUNCTION_PREFIX ## _free(MATRIX_TYPE *m)

#define MATRIX_EXP(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, EXPFUNCTION) \
    void FUNCTION_PREFIX ## _exp(MATRIX_TYPE *m) \
    { \
        size_t i, dim = m->size; \
        TYPE *M = m->M; \
\
        for(i = 0; i < dim*dim; i++) \
            M[i] = EXPFUNCTION(M[i]); \
    } 

#define MATRIX_EXP_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) void FUNCTION_PREFIX ## _exp(MATRIX_TYPE *m)

#define MATRIX_FROEBENIUS(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, SQRTFUNCTION) \
    TYPE FUNCTION_PREFIX ## _froebenius(MATRIX_TYPE *M) \
    { \
        size_t i,j, dim = M->size; \
        TYPE norm = 0; \
\
        for(i = 0; i < dim; i++) \
            for(j = 0; j < dim; j++) \
                norm += pow_2(matrix_get(M, i,j)); \
\
        return SQRTFUNCTION(norm); \
    }

#define MATRIX_FROEBENIUS_HEADER(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) TYPE FUNCTION_PREFIX ## _froebenius(MATRIX_TYPE *M)

#define MATRIX_LOGDET_LAPACK(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, ABS_FUNCTION, COPYSIGN_FUNCTION, SQRT_FUNCTION, LOG_FUNCTION) \
    TYPE FUNCTION_PREFIX ## _logdet(MATRIX_TYPE *M) \
    { \
        double det = 0; \
        int i,dim = M->size; \
        int IPIV[dim]; \
        LAPACKE_dgetrf(LAPACK_COL_MAJOR, dim, dim, M->M, dim,IPIV); \
\
        for(i = 0; i < dim; i++) \
            det += LOG_FUNCTION(ABS_FUNCTION(matrix_get(M, i, i))); \
        return det; \
    }


#define MATRIX_LOGDET(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, ABS_FUNCTION, COPYSIGN_FUNCTION, SQRT_FUNCTION, LOG_FUNCTION) \
    TYPE FUNCTION_PREFIX ## _logdet(MATRIX_TYPE *M) \
    { \
        size_t i, j, n, dim = M->size; \
        TYPE det = 0; \
\
        for(j = 0; j < dim-1; j++) \
            for(i = j+1; i < dim; i++) \
            {\
                TYPE c,s, Mij = matrix_get(M, i, j); \
\
                if(Mij != 0) \
                { \
                    TYPE a = matrix_get(M, j, j); \
                    TYPE b = Mij; \
 \
                    if(b == 0) \
                    { \
                        c = COPYSIGN_FUNCTION(1,a); \
                        s = 0; \
                    } \
                    else if(a == 0) \
                    { \
                        c = 0; \
                        s = -COPYSIGN_FUNCTION(1, b); \
                    } \
                    else if(ABS_FUNCTION(b) > ABS_FUNCTION(a)) \
                    { \
                        TYPE t = a/b; \
                        TYPE u = COPYSIGN_FUNCTION(SQRT_FUNCTION(1+t*t),b); \
                        s = -1/u; \
                        c = -s*t; \
                    } \
                    else \
                    { \
                        TYPE t = b/a; \
                        TYPE u = COPYSIGN_FUNCTION(SQRT_FUNCTION(1+t*t),a); \
                        c = 1/u; \
                        s = -c*t; \
                    } \
 \
                    for(n = 0; n < dim; n++) \
                    { \
                        TYPE Min = matrix_get(M, i, n); \
                        TYPE Mjn = matrix_get(M, j, n); \
 \
                        matrix_set(M, i, n, c*Min + s*Mjn); \
                        matrix_set(M, j, n, c*Mjn - s*Min); \
                    } \
 \
                    matrix_set(M, i, j, 0); \
                } \
            } \
 \
        for(i = 0; i < dim; i++) \
            det += LOG_FUNCTION(ABS_FUNCTION(matrix_get(M, i, i))); \
        return det; \
    }

#define MATRIX_LOGDET_HEADER(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) TYPE FUNCTION_PREFIX ## _logdet(MATRIX_TYPE *M)

#define MATRIX_ABSMAX(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, ABS_FUNCTION) \
    TYPE FUNCTION_PREFIX ## _absmax(MATRIX_TYPE *M) \
    { \
        size_t i,j, dim = M->size; \
        TYPE max = ABS_FUNCTION(matrix_get(M, 0,0)); \
 \
        for(i = 0; i < dim; i++) \
            for(j = 0; j < dim; j++) \
                max = MAX(max, ABS_FUNCTION(matrix_get(M, i,j))); \
 \
        return max; \
    }

#define MATRIX_ABSMAX_HEADER(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) TYPE FUNCTION_PREFIX ## _absmax(MATRIX_TYPE *M) \

#define MATRIX_ABSMIN(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, ABS_FUNCTION) \
    TYPE FUNCTION_PREFIX ## _absmin(MATRIX_TYPE *M) \
    { \
        size_t i,j, dim = M->size; \
        TYPE max = ABS_FUNCTION(matrix_get(M, 0,0)); \
 \
        for(i = 0; i < dim; i++) \
            for(j = 0; j < dim; j++) \
                max = MIN(max, ABS_FUNCTION(matrix_get(M, i,j))); \
 \
        return max; \
    }

#define MATRIX_ABSMIN_HEADER(FUNCTION_PREFIX, MATRIX_TYPE, TYPE) TYPE FUNCTION_PREFIX ## _absmin(MATRIX_TYPE *M) \

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
#define MATRIX_BALANCE(FUNCTION_PREFIX, MATRIX_TYPE, TYPE, ABS_FUNCTION) \
    void FUNCTION_PREFIX ## _balance(MATRIX_TYPE *A) \
    { \
        size_t i,j; \
        const size_t N = A->size; \
        TYPE D[N]; \
        int not_converged = 1; \
 \
        /* initialize D to the identity matrix */ \
        for(i = 0; i < N; i++) \
            D[i] = 1; \
 \
        while(not_converged) \
        { \
            TYPE g, f, s; \
            TYPE row_norm, col_norm; \
 \
            not_converged = 0; \
 \
            for (i = 0; i < N; ++i) \
            { \
                row_norm = 0; \
                col_norm = 0; \
 \
                for (j = 0; j < N; ++j) \
                    if (j != i) \
                    { \
                      col_norm += ABS_FUNCTION(matrix_get(A, j, i)); \
                      row_norm += ABS_FUNCTION(matrix_get(A, i, j)); \
                    } \
 \
                if ((col_norm == 0.0) || (row_norm == 0.0)) \
                  continue; \
 \
                g = row_norm / FLOAT_RADIX; \
                f = 1.0; \
                s = col_norm + row_norm; \
 \
                /* \
                 * find the integer power of the machine radix which \
                 * comes closest to balancing the matrix \
                 */ \
                while (col_norm < g) \
                { \
                    f *= FLOAT_RADIX; \
                    col_norm *= FLOAT_RADIX_SQ; \
                } \
 \
                g = row_norm * FLOAT_RADIX; \
 \
                while (col_norm > g) \
                { \
                    f /= FLOAT_RADIX; \
                    col_norm /= FLOAT_RADIX_SQ; \
                } \
 \
                if ((row_norm + col_norm) < 0.95 * s * f) \
                { \
                    int k; \
                    not_converged = 1; \
 \
                    g = 1.0 / f; \
 \
                    /* \
                     * apply similarity transformation D, where \
                     * D_{ij} = f_i * delta_{ij} \
                     */ \
 \
                    /* multiply by D^{-1} on the left */ \
                    for(k = 0; k < N; k++) \
                        matrix_set(A, i,k, g*matrix_get(A,i,k)); \
 \
 \
                    /* multiply by D on the right */ \
                    for(k = 0; k < N; k++) \
                        matrix_set(A, k,i, f*matrix_get(A,k,i)); \
 \
                    /* keep track of transformation */ \
                    for(k = 0; k < N; k++) \
                        D[k] *= f; \
                } \
            } \
        } \
    }

#define MATRIX_BALANCE_HEADER(FUNCTION_PREFIX, MATRIX_TYPE) void FUNCTION_PREFIX ## _balance(MATRIX_TYPE *A)

#define matrix_get(m, i, j)   (m->M[(i)*m->size+(j)])
#define matrix_set(m, i, j,v) (m->M[(i)*m->size+(j)]=v)

MATRIX_ALLOC_HEADER(matrix, matrix_t);
MATRIX_FREE_HEADER (matrix, matrix_t);
MATRIX_EXP_HEADER  (matrix, matrix_t);
MATRIX_FROEBENIUS_HEADER(matrix, matrix_t, double);
MATRIX_LOGDET_HEADER(matrix, matrix_t, double);
MATRIX_ABSMIN_HEADER(matrix, matrix_t, double);
MATRIX_ABSMAX_HEADER(matrix, matrix_t, double);
MATRIX_BALANCE_HEADER(matrix, matrix_t);

MATRIX_ALLOC_HEADER(matrix_char, matrix_char_t);
MATRIX_FREE_HEADER (matrix_char, matrix_char_t);

MATRIX_ALLOC_HEADER(matrix_edouble, matrix_edouble_t);
MATRIX_FREE_HEADER (matrix_edouble, matrix_edouble_t);
MATRIX_EXP_HEADER  (matrix_edouble, matrix_edouble_t);
MATRIX_FROEBENIUS_HEADER(matrix_edouble, matrix_edouble_t, edouble);
MATRIX_LOGDET_HEADER    (matrix_edouble, matrix_edouble_t, edouble);
MATRIX_ABSMIN_HEADER    (matrix_edouble, matrix_edouble_t, edouble);
MATRIX_ABSMAX_HEADER    (matrix_edouble, matrix_edouble_t, edouble);
MATRIX_BALANCE_HEADER   (matrix_edouble, matrix_edouble_t);

#endif
