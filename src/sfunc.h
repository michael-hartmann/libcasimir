#ifndef __SFUNC_H
#define __SFUNC_H

#define M_LOGPI 1.1447298858494002

// abbrevations for functions
#define lngamma(x) (lgamma(x))
#define pow_2(x) (x*x)
#define lnfac(x) (lgamma(1+x))

double logadd(double a, double b);
double logadd_s(double a, int sign_a, double b, int sign_b, int *sign);

double inline binom(int n, int k);

double bessel_lnInu(const int n, const double x);
double bessel_lnKnu(const int n, const double x);
void bessel_lnInuKnu(const int nu, const double x, double *lnInu, double *lnKnu);

#endif
