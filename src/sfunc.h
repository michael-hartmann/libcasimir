#ifndef __SFUNC_H
#define __SFUNC_H

#define M_LOGPI 1.1447298858494002

// abbrevations for functions
#define lngamma(x) (lgamma(x))
#define pow_2(x) (x*x)
#define lnfac(x) (lgamma(1+x))

#include <gsl/gsl_sf.h>

/*
#define bessel_lnInu(n,x) (log(gsl_sf_bessel_Inu(n+0.5,x)))
#define bessel_lnKnu(n,x) (gsl_sf_bessel_lnKnu(n+0.5,x))
*/
double bessel_lnInu(const int n, const double x);
double bessel_lnKnu(const int n, const double x);
double bessel_lnInuKnu(const int nu, const double x, double *lnInu, double *lnKnu);

#endif
