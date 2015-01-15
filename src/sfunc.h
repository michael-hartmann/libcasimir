#ifndef __SFUNC_H
#define __SFUNC_H

#define M_LNPI 1.1447298858494002
#define M_LN4  1.3862943611198906

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

double ln_doublefact(int n);

double plm_lnPlm(int l, int m, double x, int *sign);
double plm_Plm(int l, int m, double x);
double plm_lndPlm(int l, int m, double x, int *sign);
double plm_dPlm(int l, int m, double x);

void plm_PlmPlm(int l1, int l2, int m, double x, plm_combination_t *res);

#endif
