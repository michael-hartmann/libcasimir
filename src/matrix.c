#include <stdio.h>

#include "libcasimir.h"
#include "sfunc.h"
#include "matrix.h"

MATRIX_ALLOC(matrix_edouble, matrix_edouble_t, edouble);
MATRIX_FREE (matrix_edouble, matrix_edouble_t);
MATRIX_EXP       (matrix_edouble, matrix_edouble_t, edouble, expq);
MATRIX_FROEBENIUS(matrix_edouble, matrix_edouble_t, edouble, sqrtq);
MATRIX_LOGDET    (matrix_edouble, matrix_edouble_t, edouble, fabsq, copysignq, sqrtq, logq);
MATRIX_ABSMIN    (matrix_edouble, matrix_edouble_t, edouble, fabsq);
MATRIX_ABSMAX    (matrix_edouble, matrix_edouble_t, edouble, fabsq);
MATRIX_BALANCE   (matrix_edouble, matrix_edouble_t, edouble, fabsq);
