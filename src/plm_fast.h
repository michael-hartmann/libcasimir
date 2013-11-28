#ifndef __PLM_FAST_H
#define __PLM_FAST_H

#include <gsl/gsl_matrix.h>

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

double Yl1mYl2m(int l1, int l2, int m, double x);
double dYl1mdYl2m(int l1, int l2, int m, double x);
double Yl1mdYl2m(int l1, int l2, int m, double x);
void YYY(int l1, int l2, int m, double x, double *yl1myl2m, double *dyl1mdyl2m, double *yl1mdyl2m, double *dyl1myl2m);

#endif
