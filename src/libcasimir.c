#define _GNU_SOURCE

#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>

#include "integration.h"
#include "libcasimir.h"
#include "matrix.h"
#include "sfunc.h"
#include "utils.h"

#define FACTOR_LMAX 5
#define HBARC 3.161526510740123e-26
#define KB    1.3806488e-23

#define EPS_PRECISION 1e-16

char CASIMIR_COMPILE_INFO[255];

const char *casimir_compile_info(void)
{
    snprintf(CASIMIR_COMPILE_INFO, sizeof(CASIMIR_COMPILE_INFO)/sizeof(char), "Compiler %s, using %s", COMPILER, CASIMIR_ARITHMETICS);
    return CASIMIR_COMPILE_INFO;
}

/* casimir_lnLambda 
 * This function returns the logarithm of the Λ prefactor for given l1,l2,m.
 *
 * The values are computed using the lngamma function, so that nominator
 * and denominator of the term
 * (l1-m)!*(l2-m)!/(l1+m)!/(l2+m) 
 * don't overflow.
 *
 * Lambda(l1,l2,m) = -2*N_{l1,m}*N_{l2,m} / sqrt(l1*(l1+1)*l2*(l2+1))
 *
 * Restrictions: l1,l2,m integers, l1,l2>=1, l1,l2 >= m
 * Symmetries: Λ(l1,l2,m) = Λ(l2,l1,m)
 */
double inline casimir_lnLambda(int l1, int l2, int m, int *sign)
{
    if(sign != NULL)
        *sign = -1;
    return M_LN2 + (log(2.*l1+1)+log(2*l2+1)-log(4)-log(l1)-log(l1+1)-log(l2)-log(l2+1)+lnfac(l1-m)+lnfac(l2-m)-lnfac(l1+m)-lnfac(l2+m))/2.0;
}


/*
 * Convert the Free Energy F in SI units to Free Energy F in units of Script/(hbar*c)
 */
double casimir_F_SI_to_scaled(double F_SI, double ScriptL_SI)
{
    return ScriptL_SI/(HBARC)*F_SI;
}

/*
 * Convert the Free Energie F in units of ScriptL/(hbar*c) to Free energy F in SI units
 */
double casimir_F_scaled_to_SI(double F, double ScriptL_SI)
{
    return HBARC/ScriptL_SI*F;
}

/*
 * Convert the temperature T in Kelvin to temperature in units of 2pi*kb*ScriptL/(hbar*c)
 */
double casimir_T_SI_to_scaled(double T_SI, double ScriptL_SI)
{
    return 2*M_PI*KB*ScriptL_SI/HBARC*T_SI;
}

/*
 * Convert the temperature T in units of 2pi*kb*ScriptL/(hbar*c) to Kelvin */
double casimir_T_scaled_to_SI(double T, double ScriptL_SI)
{
    return HBARC/(2*M_PI*KB*ScriptL_SI)*T;
}

/* Function Xi
 * This function returns the Ξ prefactor for given l1,l2,m
 *
 * lgamma is used to prevent overflows - like in Λ.
 *
 * Restrictions: l1,l2,m integers, l1,l2>=1, l1,l2 >= m
 */
double casimir_lnXi(int l1, int l2, int m, int *sign)
{
    *sign = pow(-1, l2);
    return (log(2*l1+1)+log(2*l2+1)-lnfac(l1-m)-lnfac(l2-m)-lnfac(l1+m)-lnfac(l2+m)-log(l1)-log(l1+1)-log(l2)-log(l2+1))/2.0 \
           +lnfac(2*l1)+lnfac(2*l2)+lnfac(l1+l2)-M_LOG4*(2*l1+l2+1)-lnfac(l1-1)-lnfac(l2-1);
}


/*
 * The Casimir class provides methods to calculate the free energy of the
 * Casimir effect in the plane-sphere geometry for perfect reflectors.
 *
 * The main goal is to calculate the free energy and derived quantities.
 * Create a new Casimir object.
 * 
 * RbyScriptL: ratio of radius of sphere and distance between center of sphere
 *             and plate: R/(R+L)
 * T:          Temperature
 *
 * Restrictions: T > 0, 0 < RbyScriptL < 1
 */
int casimir_init(casimir_t *self, double RbyScriptL, double T)
{
    if(RbyScriptL < 0 || RbyScriptL >= 1 || T < 0)
        return -1;

    self->lmax = (int)ceil(RbyScriptL*FACTOR_LMAX);

    self->T           = T;
    self->RbyScriptL  = RbyScriptL;
    self->precision   = EPS_PRECISION;
    self->extrapolate = 0;
    self->verbose     = 0;
    self->cores       = 1;
    self->threads     = NULL;

    return 0;
}

int casimir_get_extrapolate(casimir_t *self)
{
    return self->extrapolate;
}

int casimir_set_extrapolate(casimir_t *self, int extrapolate)
{
    self->extrapolate = extrapolate ? 1 : 0;
    return 1;
}

int casimir_get_cores(casimir_t *self)
{
    return self->cores;
}

int casimir_set_cores(casimir_t *self, int cores)
{
    if(cores < 1)
        return 0;

    self->cores = cores;
    self->threads = xrealloc(self->threads, cores*sizeof(pthread_t));

    return 1;
}

/*
 * Set maximum value for l. lmax must be positive.
 */
int casimir_set_lmax(casimir_t *self, int lmax)
{
    if(lmax > 0)
        return 0;

    self->lmax = lmax;
    return 1;
}

int casimir_get_lmax(casimir_t *self)
{
    return self->lmax;
}

int casimir_get_verbose(casimir_t *self)
{
    return self->verbose;
}

int casimir_set_verbose(casimir_t *self, int verbose)
{
    self->verbose = verbose ? 1 : 0;
    return 1;
}

double casimir_get_precision(casimir_t *self)
{
    return self->precision;
}

int casimir_set_precision(casimir_t *self, double precision)
{
    if(precision <= 0)
        return 0;

    self->precision = precision;
    return 1;
}

/*
 * Free memory for casimir object
 */
void casimir_free(casimir_t *self)
{
    if(self->threads != NULL)
    {
        xfree(self->threads);
        self->threads = NULL;
    }
}


/*
 * For small x<<1 a_l will scale as
 *      a_l(x) ~ a0*(x/2)^(2l+1)
 *
 * For small x<<1 b_l will scale as
 *      b_l(x) ~ b0*(x/2)^(2l+1)
 *
 * This method returns the prefactor a0, b0. The signs of a0 and b0 are stored
 * in sign_a0 and sign_b0.
 */
void casimir_lna0_lnb0(int l, double *a0, int *sign_a0, double *b0, int *sign_b0)
{
    *sign_a0 = pow(-1, l);
    *sign_b0 = pow(-1, l+1);
    *b0 = M_LOGPI-lngamma(l+0.5)-lngamma(l+1.5);
    *a0 = *b0+log1p(1.0/l);
}

double casimir_lna(int l, const double arg, int *sign)
{
    double nominator, denominator, frac, ret;
    double lnKlp,lnKlm,lnIlm,lnIlp;
    double prefactor;
    double lnfrac = log(arg)-log(l);

    /* we could do both calculations together. but it doesn't cost much time -
     * so why bother? 
     */
    bessel_lnInuKnu(l-1, arg, &lnIlm, &lnKlm);
    bessel_lnInuKnu(l,   arg, &lnIlp, &lnKlp);

    prefactor = M_LOGPI-M_LN2+lnIlp-lnKlp;
    *sign = pow(-1, l+1);

    // nominator
    {
        frac = exp(lnfrac+lnIlm-lnIlp);
        if(frac < 1)
            nominator = log1p(fabs(-frac));
        else
        {
            if(frac > 1)
                *sign *= -1;

            nominator = log(fabs(1-frac));
        }
    }
    // denominator
    {
        frac = exp(lnfrac+lnKlm-lnKlp);
        if(frac < 1)
            denominator = log1p(frac);
        else
            denominator = log(1+frac);
    }

    ret = prefactor+nominator-denominator;

    assert(!isnan(ret));
    assert(!isinf(ret));

    return ret;
}

/*        
 * Returns the coefficient b_l for reflection on the sphere
 *
 * Restrictions: l integer, l>=1, xi>0
 */        
double casimir_lnb(int l, const double arg, int *sign)
{
    double lnInu, lnKnu, ret;

    bessel_lnInuKnu(l, arg, &lnInu, &lnKnu);
    *sign = pow(-1, l+1);

    ret = M_LOGPI-M_LN2+lnInu-lnKnu;

    assert(!isnan(ret));
    assert(!isinf(ret));

    return ret;
}

/*
 * Initialize the mie cache.
 * This function must be called before any call to casimir_mie_cache_alloc
 */
void casimir_mie_cache_init(casimir_mie_cache_t *cache, double arg)
{
    cache->al = cache->bl = NULL;
    cache->al_sign = cache->bl_sign = NULL;
    cache->lmax = 0;
    cache->arg = arg;
}

/*
 * Allocate memory for the Mie-coefficients a_l and b_l
 */
int casimir_mie_cache_alloc(casimir_t *self, casimir_mie_cache_t *cache, int lmax)
{
    int l;
    double arg = cache->arg;

    if(arg == 0)
    {
        cache->al = cache->bl = NULL;
        cache->al_sign = cache->bl_sign = NULL;
        return 0;
    }

    cache->al      = (double *)xrealloc(cache->al,      (lmax+1)*sizeof(double));
    cache->bl      = (double *)xrealloc(cache->bl,      (lmax+1)*sizeof(double));
    cache->al_sign =    (int *)xrealloc(cache->al_sign, (lmax+1)*sizeof(int));
    cache->bl_sign =    (int *)xrealloc(cache->bl_sign, (lmax+1)*sizeof(int));

    cache->al[0] = cache->bl[0] = 0;
    for(l = MAX(1,cache->lmax); l <= lmax; l++)
    {
        cache->al[l] = casimir_lna(l,arg, &cache->al_sign[l]);
        cache->bl[l] = casimir_lnb(l,arg, &cache->bl_sign[l]);
    }
    cache->lmax = lmax;

    return 1;
}

/*
 * Free memory of cache.
 */
void casimir_mie_cache_free(casimir_mie_cache_t *cache)
{
    if(cache->al != NULL)
        xfree(cache->al);
    if(cache->bl != NULL)
        xfree(cache->bl);
    if(cache->al_sign != NULL)
        xfree(cache->al_sign);
    if(cache->bl_sign != NULL)
        xfree(cache->bl_sign);

    cache->al = cache->bl = NULL;
}

/* Sum len numbers in value.
   The idea is: To avoid a loss of significance, we sum beginning with smallest
   number and add up in increasing order
*/
static double _sum(double values[], size_t len)
{
    int i;
    double sum = 0;

    for(i = len-1; i > 0; i--)
        sum += values[i];

    sum += values[0]/2;

    return sum;
}

static void *_thread_f(void *p)
{
    casimir_thread_t *r = (casimir_thread_t *)p;
    r->value = casimir_F_n(r->self, r->n, &r->nmax);
    return r;
}

static pthread_t *_start_thread(casimir_thread_t *r)
{
    pthread_t *t = xmalloc(sizeof(pthread_t));
    pthread_create(t, NULL, _thread_f, (void *)r);

    return t;
}

static int _join_threads(casimir_t *self, double values[], int *ncalc)
{
    int i, joined = 0, running = 0;
    casimir_thread_t *r;
    pthread_t **threads = self->threads;

    for(i = 0; i < self->cores; i++)
    {
        if(threads[i] != NULL)
        {
            running++;

            if(pthread_tryjoin_np(*threads[i], (void *)&r) == 0)
            {
                joined++;

                if(r->n > *ncalc)
                    *ncalc = r->n;

                values[r->n] = r->value;
                xfree(r);
                xfree(threads[i]);
                threads[i] = NULL;

                if(self->verbose)
                    fprintf(stderr, "# n=%d, value=%.15g\n", r->n, values[r->n]);
            }
        }
    }

    if(running == 0)
        return -1;

    return joined;
}

double casimir_F_n(casimir_t *self, const int n, int *mmax)
{
    double precision = self->precision;
    double TRbyScriptL = self->T*self->RbyScriptL;
    casimir_mie_cache_t cache;
    double sum_n = 0;
    int m;
    const int lmax = self->lmax;
    double values[lmax+1];

    for(m = 0; m <= lmax; m++)
        values[m] = 0;

    casimir_mie_cache_init(&cache, n*TRbyScriptL);
    casimir_mie_cache_alloc(self, &cache, self->lmax);

    for(m = 0; m <= self->lmax; m++)
    {
        values[m] = casimir_logdetD(self,n,m,&cache);

        if(self->verbose)
            fprintf(stderr, "# n=%d, m=%d, value=%.15g\n", n, m, values[m]);

        /* If F is !=0 and value/F < 1e-16, then F+value = F. The addition
         * has no effect.
         * As for larger m value will be even smaller, we can skip the
         * summation here. 
         */
        sum_n = _sum(values, lmax+1);
        if(values[0] != 0 && fabs(values[m]/sum_n) < precision)
            break;
    }

    casimir_mie_cache_free(&cache);

    if(self->verbose)
        fprintf(stderr, "# n=%d, value=%.15g\n", n, sum_n);

    if(mmax != NULL)
        *mmax = m;

    return sum_n;
}

/*
 * Calculate free energy.

 * Restrictions: nmax integer, nmax >= 0
 */
double casimir_F(casimir_t *self, int *nmax)
{
    int i, n = 0;
    double sum_n = 0;
    const double precision = self->precision;
    double *values = NULL;
    size_t len = 0;
    int ncalc = 0;
    const int cores = self->cores;
    pthread_t **threads = self->threads;

    if(cores > 1)
        for(i = 0; i < cores; i++)
            threads[i] = NULL;

    /* So, here we sum up all m and n that contribute to F.
     * So, what do we do here?
     *
     * We want to evaluate
     *      \sum_{n=0}^\infty \sum{m=0}^{l_max} log(det(1-M)),
     * where the terms for n=0 and m=0 are weighted with a factor 1/2.
     */
    while(1)
    {
        if(n >= len)
        {
            const int delta = MAX(512, cores);

            values = (double *)xrealloc(values, (len+delta)*sizeof(double));

            for(i = len; i < len+delta; i++)
                values[i] = 0;

            len += delta;
        }

        if(cores > 1)
        {
            casimir_thread_t *r;

            for(i = 0; i < cores; i++)
            {
                if(threads[i] == NULL)
                {
                    r = (casimir_thread_t *)xmalloc(sizeof(casimir_thread_t));

                    r->self  = self;
                    r->n     = n++;
                    r->value = 0;
                    r->nmax  = 0;

                    threads[i] = _start_thread(r);
                }
            }

            if(_join_threads(self, values, &ncalc) == 0)
                usleep(CASIMIR_IDLE);
        }
        else
        {
            values[n] = casimir_F_n(self, n, NULL);

            ncalc = n;
            n++;
        }

        if(values[0] != 0)
        {
            if(fabs(values[ncalc]/(2*values[0])) < precision)
            {
                if(cores > 1)
                    while(_join_threads(self, values, &ncalc) != -1)
                        usleep(CASIMIR_IDLE);

                sum_n = _sum(values, n);
                if(self->extrapolate && n > 20)
                {
                    double r1 = values[n-1]/values[n-2];
                    double r2 = values[n-2]/values[n-3];
                    double r3 = values[n-3]/values[n-4];
                    double r4 = values[n-4]/values[n-5];
                    double r5 = values[n-5]/values[n-6];
                    double r  = (r1+r2+r3+r4+r5)/5;
                    sum_n += values[n-1]*r/(1-r);
                }
                /* get out of here */
                if(nmax != NULL)
                    *nmax = n-1; // we calculated n term from n=0,...,nmax=n-1

                if(values != NULL)
                    xfree(values);

                return self->T/M_PI*sum_n;
            }
        }
    }
}

double casimir_logdetD0(casimir_t *self, int m, double *logdet_EE, double *logdet_MM)
{
    int l1,l2,min,max,dim;
    double lnRbyScriptL = log(self->RbyScriptL);
    double value_EE, value_MM;

    min = MAX(m,1);
    max = self->lmax;

    dim = (max-min+1);

    matrix_edouble_t *EE = matrix_edouble_alloc(dim);
    matrix_edouble_t *MM = matrix_edouble_alloc(dim);

    /* calculate the logarithm of the matrix elements of D */
    for(l1 = min; l1 <= max; l1++)
        for(l2 = min; l2 <= max; l2++)
        {
            /* i: row of matrix, j: column of matrix */
            const int i = l1-min, j = l2-min;
            int sign_a0, sign_b0, sign_xi;
            double lna0, lnb0;
            double lnXiRL = casimir_lnXi(l1,l2,m,&sign_xi)+(2*l1+1)*lnRbyScriptL;
            casimir_lna0_lnb0(l1, &lna0, &sign_a0, &lnb0, &sign_b0);

            matrix_set(EE, i,j, (l1 == l2 ? 1 : 0) - sign_xi*sign_a0*expq(lna0+lnXiRL));
            matrix_set(MM, i,j, (l1 == l2 ? 1 : 0) - sign_xi*sign_a0*expq(lnb0+lnXiRL));
        }

    /* balance the matrix */
    matrix_edouble_balance(EE);
    matrix_edouble_balance(MM);

    value_EE = matrix_edouble_logdet(EE);
    value_MM = matrix_edouble_logdet(MM);

    /* free space for matrices */
    matrix_edouble_free(EE);
    matrix_edouble_free(MM);

    if(logdet_EE != NULL)
        *logdet_EE = value_EE;
    if(logdet_MM != NULL)
        *logdet_MM = value_MM;

    return value_EE+value_MM;
}

/*
 * Calculate logarithm of determinant of D=1-M for given m,ξ
 *
 * Restrictions: m integer, m>=0, ξ>= 0
 */
double casimir_logdetD(casimir_t *self, int n, int m, casimir_mie_cache_t *cache)
{
    int min,max,dim,l1,l2;
    double logdet = 0;
    double nTRbyScriptL = n*self->T*self->RbyScriptL;

    min = MAX(m,1);
    max = self->lmax;

    dim = (max-min+1);

    if(n == 0)
        return casimir_logdetD0(self, m, NULL, NULL);

    matrix_edouble_t *M = matrix_edouble_alloc(2*dim);

    /* M_EE, -M_EM
       M_ME,  M_MM */
    for(l1 = min; l1 <= max; l1++)
    {
        for(l2 = min; l2 <= l1; l2++)
        {
            const int i = l1-min, j = l2-min;
            casimir_integrals_t cint;
            double lnal1 = cache->al[l1];
            double lnbl1 = cache->bl[l1];
            double lnal2 = cache->al[l2];
            double lnbl2 = cache->bl[l2];

            double al1_sign = cache->al_sign[l1];
            double bl1_sign = cache->bl_sign[l1];
            double al2_sign = cache->al_sign[l2];
            double bl2_sign = cache->bl_sign[l2];

            if(nTRbyScriptL < 1)
            {
                double lognTRbyScriptL = log(nTRbyScriptL);
                lnal1 -= (l1-l2)*lognTRbyScriptL;
                lnbl1 -= (l1-l2)*lognTRbyScriptL;

                lnal2 -= (l2-l1)*lognTRbyScriptL;
                lnbl2 -= (l2-l1)*lognTRbyScriptL;
            }

            casimir_integrate(&cint, l1, l2, m, n*self->T);

            assert(!isnan(cint.logA));
            assert(!isnan(cint.logB));
            assert(!isnan(cint.logC));
            assert(!isnan(cint.logD));

            assert(!isinf(cint.logB));
            if(m != 0)
            {
                assert(!isinf(cint.logA));
                assert(!isinf(cint.logC));
                assert(!isinf(cint.logD));
            }

            /* EE */
            matrix_set(M, i,j, (l1 == l2 ? 1 : 0) - al1_sign*( cint.signB*expq(lnal1+cint.logB) - cint.signA*expq(lnal1+cint.logA) ));
            matrix_set(M, j,i, (l1 == l2 ? 1 : 0) - pow(-1, l1+l2)*al2_sign*( cint.signB*expq(lnal2+cint.logB) - cint.signA*expq(lnal2+cint.logA) ));

            /* MM */
            matrix_set(M, i+dim,j+dim, (l1 == l2 ? 1 : 0) - bl1_sign*( cint.signA*expq(lnbl1+cint.logA) - cint.signB*expq(lnbl1+cint.logB) ));
            matrix_set(M, j+dim,i+dim, (l1 == l2 ? 1 : 0) - pow(-1, l1+l2)*bl2_sign*( cint.signA*expq(lnbl2+cint.logA) - cint.signB*expq(lnbl2+cint.logB) ));


            if(m != 0)
            {
                /* M_EM */
                matrix_set(M, dim+i,j, -al1_sign*( cint.signD*expq(lnal1+cint.logD) - cint.signC*expq(lnal1+cint.logC) ));
                matrix_set(M, dim+j,i, -al2_sign*pow(-1, l1+l2+1)*( cint.signC*expq(lnal2+cint.logC) - cint.signD*expq(lnal2+cint.logD) ));


                /* M_ME */
                matrix_set(M, i,dim+j, -bl1_sign*( cint.signC*expq(lnbl1+cint.logC) - cint.signD*expq(lnbl1+cint.logD) ));
                matrix_set(M, j,dim+i, -bl2_sign*pow(-1, l1+l2+1)*( cint.signD*expq(lnbl2+cint.logD) - cint.signC*expq(lnbl2+cint.logC) ));
            }
        }
    }

    if(m == 0)
    {
        size_t i,j;
        matrix_edouble_t *EE = matrix_edouble_alloc(dim);
        matrix_edouble_t *MM = matrix_edouble_alloc(dim);

        for(i = 0; i < dim; i++)
            for(j = 0; j < dim; j++)
            {
                matrix_set(EE, i,j, matrix_get(M, i,j));
                matrix_set(MM, i,j, matrix_get(M, dim+i,dim+j));
            }

        matrix_edouble_balance(MM);
        matrix_edouble_balance(EE);

        logdet = matrix_edouble_logdet(EE)+matrix_edouble_logdet(MM);

        matrix_edouble_free(EE);
        matrix_edouble_free(MM);
    }
    else
    {
        matrix_edouble_balance(M);
        logdet = matrix_edouble_logdet(M);
    }

    matrix_edouble_free(M);

    assert(!isinf(logdet));
    return logdet;
}
