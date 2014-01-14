#ifndef __CASIMIR_H
#define __CASIMIR_H

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>

#define TE 0
#define TM 1

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

#define SLA(l,n,arg) (gsl_sf_bessel_Inu((l)+0.5, (n)*(arg)) * ((l)*gsl_sf_bessel_Inu((l)+0.5, (arg))     -     (arg)*gsl_sf_bessel_Inu((l)-0.5, (arg))))
#define SLB(l,b,arg) (gsl_sf_bessel_Inu((l)+0.5, (arg))     * ((l)*gsl_sf_bessel_Inu((l)+0.5, (n)*(arg)) - (n)*(arg)*gsl_sf_bessel_Inu((l)-0.5, (arg))))
#define SLC(l,b,arg) (gsl_sf_bessel_Inu((l)+0.5, (n)*(arg)) * ((l)*gsl_sf_bessel_Knu((l)+0.5, (arg))     +     (arg)*gsl_sf_bessel_Knu((l)-0.5, (arg))))
#define SLD(l,b,arg) (gsl_sf_bessel_Knu((l)+0.5, (arg))     * ((l)*gsl_sf_bessel_Knu((l)+0.5, (n)*(arg)) - (n)*(arg)*gsl_sf_bessel_Knu((l)-0.5, (arg))))

#define lngamma(x) (gsl_sf_lngamma(x))
#define lnfac(x) (gsl_sf_lngamma(1+x))
#define IS_FINITE(x) ((!isnan(x) && !isinf(x)))

typedef struct
{
    double RbyScriptL; // R/(R+L)
    double T;
    int lmax;
    int verbose;

    double eps;
} casimir_t;

typedef struct
{
    double *al;
    double *bl;
    int lmax;
    double arg;
} casimir_mie_cache_t;

typedef struct
{
    int l1, l2, m;
    casimir_t *self;
    double nT;
    double scale;
    double(*rp)(casimir_t *self, double x, double xi);
} casimir_int_t;

typedef struct
{
    double A,B,C,D;
} casimir_integrals_t;

double casimir_lnLambda(int l1, int l2, int m);
double casimir_Lambda(int l1, int l2, int m);
double casimir_lnXi(int l1, int l2, int m, int *sign);
double casimir_Xi(int l1,int l2, int m);

double casimir_F_SI_to_scaled(double F_SI, double ScriptL_SI);
double casimir_F_scaled_to_SI(double F, double ScriptL_SI);
double casimir_T_SI_to_scaled(double T_SI, double ScriptL_SI);
double casimir_T_scaled_to_SI(double T, double ScriptL_SI);

void casimir_lna0_lnb0(int l, double *a0, int *sign_a0, double *b0, int *sign_b0);

int casimir_init(casimir_t *self, double RbyScriptL, double T);

double casimir_get_lmax(casimir_t *self);
void casimir_set_lmax(casimir_t *self, int lmax);
void casimir_set_limits(casimir_t *self, int limits);
void casimir_set_eps(casimir_t *self, double eps);
void casimir_set_verbose(casimir_t *self, int verbose);
void casimir_free(casimir_t *self);

double casimir_a(casimir_t *self, int l, double arg);
double casimir_b(casimir_t *self, int l, double arg);

double casimir_epsilon(casimir_t *self, double xi);
double casimir_rTE(casimir_t *self, double x, double xi);
double casimir_rTM(casimir_t *self, double x, double xi);

void casimir_integrands_vec(double x, void *params, double *vec, int len);

double casimir_F(casimir_t *self, int *nmax);

void casimir_mie_cache_init(casimir_mie_cache_t *cache, double arg);
int casimir_mie_cache_alloc(casimir_t *self, casimir_mie_cache_t *cache, int lmax);
void casimir_mie_cache_free(casimir_mie_cache_t *cache);

double casimir_logdetD(casimir_t *self, int n, int m, casimir_mie_cache_t *cache);
int casimir_logdet1m(gsl_matrix *M, double *logdet, int n, int m, const char *desc);

double casimir_logdetD_approx(casimir_t *self, int n, int m, casimir_mie_cache_t *cache);
#endif
