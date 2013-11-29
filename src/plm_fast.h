#ifndef __PLM_FAST_H
#define __PLM_FAST_H

#include <gsl/gsl_matrix.h>

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

#define Nlm(l,m) ( sqrt((2.0*l+1)/2. * exp(gsl_sf_lngamma(1+l-m)-gsl_sf_lngamma(1+l+m))) )
#define plm_dYlm(l,m,x) (plm_dPlm(l,m,x)*Nlm(l,m))
#define plm_Ylm(l,m,x)  (plm_Plm(l,m,x)*Nlm(l,m))

double plm_Plm(int l, int m, double x);
double plm_dPlm(int l, int m, double x);
void plm_Yl12md(int l1, int l2, int m, double x, double *Yl1m, double *Yl2m, double *dYl1m, double *dYl2m);
double plm_Yl1mYl2m(int l1, int l2, int m, double x);
double plm_dYl1mdYl2m(int l1, int l2, int m, double x);
double plm_Yl1mdYl2m(int l1, int l2, int m, double x);

#endif
