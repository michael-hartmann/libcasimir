#ifndef __SFUNC_H
#define __SFUNC_H

#define M_LOGPI 1.1447298858494002
#define M_LOG4  1.3862943611198906

// abbrevations for functions
#define lngamma(x) (lgamma(x))
#define pow_2(x) (x*x)
#define lnfac(x) (lgamma(1+x))

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) ((((a))>((b)))?((a)):((b)))
#endif

double inline logadd(const double a, const double b);
double inline logadd_m(const double list[], size_t len);
double inline logadd_s(const double a, const int sign_a, const double b, const int sign_b, int *sign);
double inline logadd_ms(const double list[], const char signs[], const size_t len, int *sign);

double inline lbinom(int n, int k);
double inline binom(int n, int k);

double bessel_lnInu(const int n, const double x);
double bessel_lnKnu(const int n, const double x);
void bessel_lnInuKnu(const int nu, const double x, double *lnInu, double *lnKnu);

int round2up(int x);

double linspace(double start, double stop, int N, int i);
double logspace(double start, double stop, int N, int i);

#endif
