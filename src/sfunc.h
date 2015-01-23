#ifndef __SFUNC_H
#define __SFUNC_H

#include "edouble.h"

#define M_LNPI 1.1447298858494002
#define M_LN4  1.3862943611198906

// abbrevations for functions
#define lngamma(x) (lgamma(x))
#define pow_2(x) (x*x)
#define lnfac(x)  (gamma(1+x))
#define lnfacq(x) (gammaq(1+x))

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) ((((a))>((b)))?((a)):((b)))
#endif

typedef struct {
    double lnPl1mPl2m;
    int sign_Pl1mPl2m;

    double lndPl1mPl2m;
    int sign_dPl1mPl2m;

    double lnPl1mdPl2m;
    int sign_Pl1mdPl2m;

    double lndPl1mdPl2m;
    int sign_dPl1mdPl2m;
} plm_combination_t;

double inline logadd_m(const double list[], size_t len);
edouble inline logadd_s(const edouble a, const int sign_a, const edouble b, const int sign_b, int *sign);

double inline logadd_ms(const double list[], const int signs[], const size_t len, int *sign);
edouble inline logadd_msq(const edouble list[], const int signs[], const size_t len, int *sign);

double inline lbinom(int n, int k);
double inline binom(int n, int k);

edouble bessel_lnInu(const int n, const edouble x);
edouble bessel_lnKnu(const int n, const edouble x);
void bessel_lnInuKnu(int nu, const edouble x, edouble *lnInu_p, edouble *lnKnu_p);

double linspace(double start, double stop, int N, int i);
double logspace(double start, double stop, int N, int i);

edouble ln_doublefactq(int n);

edouble plm_lnPlm (int l, int m, edouble x, int *sign);
edouble plm_Plm   (int l, int m, edouble x);
edouble plm_lndPlm(int l, int m, edouble x, int *sign);
edouble plm_dPlm  (int l, int m, edouble x);

void plm_PlmPlm(int l1, int l2, int m, edouble x, plm_combination_t *res);

#endif
