#ifndef __LIBCASIMIR_H
#define __LIBCASIMIR_H

#include <pthread.h>

#define CASIMIR_IDLE 1000

const char *casimir_compile_info(void);

typedef struct
{
    double RbyScriptL; // R/(R+L)
    double T;
    int lmax;
    int verbose;
    int cores;
    double precision;
    pthread_t **threads;
} casimir_t;

typedef struct
{
    casimir_t *self;
    int n, nmax;
    double value;
} casimir_thread_t;

typedef struct
{
    double *al;
    int *al_sign;
    double *bl;
    int *bl_sign;
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
    double logA,logB,logC,logD;
    int signA, signB, signC, signD;
} casimir_integrals_t;

double casimir_lnLambda(int l1, int l2, int m);
double casimir_lnXi(int l1, int l2, int m, int *sign);

double casimir_F_SI_to_scaled(double F_SI, double ScriptL_SI);
double casimir_F_scaled_to_SI(double F, double ScriptL_SI);
double casimir_T_SI_to_scaled(double T_SI, double ScriptL_SI);
double casimir_T_scaled_to_SI(double T, double ScriptL_SI);

void casimir_lna0_lnb0(int l, double *a0, int *sign_a0, double *b0, int *sign_b0);

int casimir_init(casimir_t *self, double RbyScriptL, double T);

double casimir_get_lmax(casimir_t *self);
int  casimir_set_cores(casimir_t *self, int cores);
void casimir_set_lmax(casimir_t *self, int lmax);
void casimir_set_limits(casimir_t *self, int limits);
void casimir_set_precision(casimir_t *self, double precision);
void casimir_set_verbose(casimir_t *self, int verbose);
void casimir_free(casimir_t *self);

double casimir_lna(int l, const double arg, int *sign);
double casimir_lnb(int l, const double arg, int *sign);

double casimir_F_n(casimir_t *self, const int n, int *mmax);
double casimir_F(casimir_t *self, int *nmax);

void casimir_mie_cache_init(casimir_mie_cache_t *cache, double arg);
int casimir_mie_cache_alloc(casimir_t *self, casimir_mie_cache_t *cache, int lmax);
void casimir_mie_cache_free(casimir_mie_cache_t *cache);

double casimir_logdetD0(casimir_t *self, int m, double *EE, double *MM);
double casimir_logdetD(casimir_t *self, int n, int m, casimir_mie_cache_t *cache);

#endif
