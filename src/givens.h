#ifndef __GIVENS_H
#define __GIVENS_H
#include <stdio.h>

#define FLOAT_RADIX       2.0
#define FLOAT_RADIX_SQ    (FLOAT_RADIX * FLOAT_RADIX)
#define LOG_FLOAT_RADIX   M_LN2
#define LOG_FLOAT_RADIX_SQ 2*M_LN2
#define LOG_095 -0.05129329438755058

#ifdef QUAD_PRECISION
    #include <quadmath.h>

    #define __float    __float128
    #define __abs      fabsq
    #define __sqrt     sqrtq
    #define __log      logq
    #define __exp      expq
    #define __copysign copysignq
#else
    #include <math.h>

    #define __float double
    #define __abs fabs
    #define __sqrt sqrt
    #define __log log
    #define __exp exp
    #define __copysign copysign
#endif

typedef struct {
    size_t size;
    __float *M;
} matrix_t;

matrix_t *matrix_alloc(size_t size);
void matrix_free(matrix_t *m);
void matrix_fprintf(const matrix_t *m, FILE *stream, const char *format, const char *lines, const char *rows);

#define matrix_get(m, i, j)   (m->M[(i)*m->size+(j)])
#define matrix_set(m, i, j,v) (m->M[(i)*m->size+(j)]=v)

__float matrix_logdet(matrix_t *M);

void matrix_balance(matrix_t *A);
void matrix_log_balance(matrix_t *A);

void matrix_info(FILE *stream, matrix_t *M);
__float matrix_froebenius(matrix_t *M);
__float matrix_absmax(matrix_t *M);
__float matrix_absmin(matrix_t *M);

void matrix_exp(matrix_t *M);

#endif
