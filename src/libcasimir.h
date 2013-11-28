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


double gsl_matrix_trace(gsl_matrix *A);

#define SLA(l,n,arg) (gsl_sf_bessel_Inu((l)+0.5, (n)*(arg)) * ((l)*gsl_sf_bessel_Inu((l)+0.5, (arg))     -     (arg)*gsl_sf_bessel_Inu((l)-0.5, (arg))))
#define SLB(l,b,arg) (gsl_sf_bessel_Inu((l)+0.5, (arg))     * ((l)*gsl_sf_bessel_Inu((l)+0.5, (n)*(arg)) - (n)*(arg)*gsl_sf_bessel_Inu((l)-0.5, (arg))))
#define SLC(l,b,arg) (gsl_sf_bessel_Inu((l)+0.5, (n)*(arg)) * ((l)*gsl_sf_bessel_Knu((l)+0.5, (arg))     +     (arg)*gsl_sf_bessel_Knu((l)-0.5, (arg))))
#define SLD(l,b,arg) (gsl_sf_bessel_Knu((l)+0.5, (arg))     * ((l)*gsl_sf_bessel_Knu((l)+0.5, (n)*(arg)) - (n)*(arg)*gsl_sf_bessel_Knu((l)-0.5, (arg))))

typedef struct
{
    double S;
    double R;
    double L;
    double T;
    double gamma;
    double omegap;
    double hbar;
    double c;
    double kBT;
    double kBT_unscaled;
    int lmax;
    int verbose;

    double epsrel;
    double eps_n;
    int int_limits;
    gsl_integration_workspace *int_workspace;
} casimir_t;

typedef struct
{
    double *al;
    double *bl;
} casimir_mie_cache_t;

typedef struct
{
    int l1, l2, m;
    casimir_t *self;
    double xi_t;
    double(*rp)(casimir_t *self, double x, double xi);
} casimir_int_t;

double Lambda(int l1, int l2, int m);
double Xi(int l1,int l2, int m);

double a0(int l);
double b0(int l);

int casimir_init(casimir_t *self, double R, double L, double T, double omegap, double gamma);
int casimir_init_perfect(casimir_t *self, double R, double L, double T);
int casimir_init_perfect_scaled(casimir_t *self, double S, double R, double L, double T);
int casimir_init_scaled(casimir_t *self, double S, double R, double L, double T, double omegap, double gamma);

void casimir_set_lmax(casimir_t *self, int lmax);
void casimir_set_limits(casimir_t *self, int limits);
void casimir_set_eps_n(casimir_t *self, double eps_n);
void casimir_set_epsrel(casimir_t *self, double epsrel);
void casimir_free(casimir_t *self);

double casimir_a(casimir_t *self, int l, double xi);
double casimir_b(casimir_t *self, int  l, double xi);

double epsilon(casimir_t *self, double xi);
double r_TE(casimir_t *self, double x, double xi);
double r_TM(casimir_t *self, double x, double xi);

double casimir_integrate(casimir_t *self, double(callback(double,void*)), int l1, int l2, int m, int p, double xi_t);
double integrandA(double x, void *params);
double integrandB(double x, void *params);
double integrandC(double x, void *params);

double casimir_intA(casimir_t *self, int l1, int l2, int m, int p, double xi);
double casimir_intB(casimir_t *self, int l1, int l2, int m, int p, double xi);
double casimir_intC(casimir_t *self, int l1, int l2, int m, int p, double xi);
double casimir_intD(casimir_t *self, int l1, int l2, int m, int p, double xi);

double xi_n(casimir_t *self, int n);

double casimir_M_EE(casimir_t *self, int l1, int l2, int m, double xi);
double casimir_M_MM(casimir_t *self, int l1, int l2, int m, double xi);
double casimir_M_EM(casimir_t *self, int l1, int l2, int m, double xi);
double casimir_M_ME(casimir_t *self, int l1, int l2, int m, double xi);

double casimir_F(casimir_t *self, int *nmax);

int casimir_mie_cache_alloc(casimir_t *self, casimir_mie_cache_t *cache, double xi);
void casimir_mie_cache_free(casimir_mie_cache_t *cache);

double casimir_logdetD(casimir_t *self, int m, double xi, casimir_mie_cache_t *cache);
int casimir_logdet1m(gsl_matrix *M, double *logdet);

#endif
