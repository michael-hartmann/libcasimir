#include <stdio.h>
#include <stdlib.h>

#include "libcasimir.h"
#include "sfunc.h"
#include "matrix.h"

MATRIX_ALLOC(matrix, matrix_t, double);
MATRIX_FREE (matrix, matrix_t);
MATRIX_EXP  (matrix, matrix_t, double, exp);
MATRIX_FROEBENIUS (matrix, matrix_t, double, sqrt);
MATRIX_LOGDET(matrix, matrix_t, double, fabs, copysign, sqrt, log);
MATRIX_ABSMIN     (matrix, matrix_t, double, fabs);
MATRIX_ABSMAX     (matrix, matrix_t, double, fabs);
MATRIX_BALANCE    (matrix, matrix_t, double, fabs);
MATRIX_LOG_BALANCE(matrix, matrix_t, double, log);

MATRIX_ALLOC(matrix_char, matrix_char_t, char);
MATRIX_FREE (matrix_char, matrix_char_t);

#ifdef MATRIX_QUAD
    MATRIX_ALLOC(matrix_quad, matrix_quad_t, __float128);
    MATRIX_FREE (matrix_quad, matrix_quad_t);
    MATRIX_EXP       (matrix_quad, matrix_quad_t, __float128, expq);
    MATRIX_FROEBENIUS(matrix_quad, matrix_quad_t, __float128, sqrtq);
    MATRIX_LOGDET    (matrix_quad, matrix_quad_t, __float128, fabsq, copysignq, sqrtq, logq);
    MATRIX_ABSMIN    (matrix_quad, matrix_quad_t, __float128, fabsq);
    MATRIX_ABSMAX    (matrix_quad, matrix_quad_t, __float128, fabsq);
    MATRIX_BALANCE   (matrix_quad, matrix_quad_t, __float128, fabsq);
#endif
