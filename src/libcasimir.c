#define _ISOC99_SOURCE
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

#include "la.h"
#include "logdet1m.h"
//#include "plm_fast.h"
#include "integration.h"
#include "libcasimir.h"
#include "gausslaguerre.h"

#define FACTOR_LMAX 5
#define HBARC 3.161526510740123e-26
#define KB   1.3806488e-23

#define EPS_PRECISION 1e-10

void handler(const char *reason, const char *file, int line, int gsl_errno);

/*
 * This is the error handler for the Gnu scientifiy library.
 * In case there are any errors, print a debug message and exit
 */
void handler(const char *reason, const char *file, int line, int gsl_errno)
{
    fprintf(stderr, "# %s in %s:%d (gsl_errno: %d)\n", reason, file, line, gsl_errno);
    exit(1);
}

/* casimir_lnLambda 
 * This function returns the logarithm of the Λ prefactor for given l1,l2,m.
 *
 * The values are computed using the lngamma function, so that nominator
 * and denominator of the term
 * (l1-m)!*(l2-m)!/(l1+m)!/(l2+m) 
 * don't overflow.
 *
 * Lambda(l1,l2,m) = N_{l1,m}*N_{l2,m} / sqrt(4*l1*(l1+1)*l2*(l2+1))
 *
 * Restrictions: l1,l2,m integers, l1,l2>=1, l1,l2 >= m
 * Symmetries: Λ(l1,l2,m) = Λ(l2,l1,m)
 */
double inline casimir_lnLambda(int l1, int l2, int m)
{
    return (log((2.*l1+1)*(2*l2+1)/(4*l1*(l1+1)*l2*(l2+1)))+gsl_sf_lngamma(1+l1-m)+gsl_sf_lngamma(1+l2-m)-gsl_sf_lngamma(1+l1+m)-gsl_sf_lngamma(1+l2+m))/2.0;
}

/* casimir_Lambda
 * This function returns the Λ prefactor for given l1,l2,m.
 *
 * See casimir_lnLambda
 */
double inline casimir_Lambda(int l1, int l2, int m)
{
    return exp(casimir_lnLambda(l1, l2, m));
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

/* casimir_logdet1m
 * This function returns returns log(det(1-M))
 *
 * The value is calculated either using Mercator series, if |M| < 1 or
 * by calculating the eigenvalues of M
 */
int casimir_logdet1m(gsl_matrix *M, double *logdet, int n, int m, const char *desc)
{
    int ret;
    // if |M| < 1 use Mercator series, otherwise use eigenvalues
    /*
    if(norm <= 0.97)
        ret = logdet1m_taylor(M, logdet);
    else
        ret = logdet1m_eigenvalues(M, logdet);
    */

    //la_matrix_info(stderr, M);
    ret = logdet1m_eigenvalues(M, logdet);
    if(ret != 0)
    {
        int i,j;
        int nans = 0, infs = 0, zeros = 0;
        double norm = la_norm_froebenius(M);
        for(i = 0; i < M->size1; i++)
            for(j = 0; j < M->size1; j++)
            {
                double elem = gsl_matrix_get(M, i,j);
                if(isinf(elem)) infs++;
                if(isnan(elem)) nans++;
                if(elem == 0)   zeros++;
            }

        fprintf(stderr, "# Can't calculate log(det(1-%s)): %d (n=%d, m=%d, norm=%g, %d zeros, %d nans, %d infs, dim %dx%d)\n", \
                        desc, ret, n, m, norm, zeros, nans, infs, (int)M->size1, (int)M->size2);
    }

    return ret;
}

/* Function Xi
 * This function returns the Ξ prefactor for given l1,l2,m
 *
 * lgamma is used to prevent overflows - like in Λ.
 *
 * Restrictions: l1,l2,m integers, l1,l2>=1, l1,l2 >= m
 */
double casimir_Xi(int l1, int l2, int m)
{
    return pow(-1, l2) * exp( \
        (log(2*l1+1)+log(2*l2+1)-lnfac(l1-m)-lnfac(l2-m)-lnfac(l1+m)-lnfac(l2+m)-log(l1)-log(l1+1)-log(l2)-log(l2+1))/2 \
        +lnfac(2*l1)+lnfac(2*l2)+lnfac(l1+l2)-log(4)*(2*l1+l2+1)-lnfac(l1-1)-lnfac(l2-1)
    );
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

    gsl_set_error_handler(&handler);

    self->lmax = (int)ceil(RbyScriptL*FACTOR_LMAX);

    self->T          = T;
    self->RbyScriptL = RbyScriptL;
    self->eps_n      = 1e-20;
    self->verbose    = 0;

    return 0;
}

/*
 * Set maximum value for l
 */
void casimir_set_lmax(casimir_t *self, int lmax)
{
    self->lmax = lmax;
}

void casimir_set_verbose(casimir_t *self, int verbose)
{
    self->verbose = verbose;
}

void casimir_set_eps_n(casimir_t *self, double eps_n)
{
    self->eps_n = eps_n;
}

/*
 * Free memory for casimir object
 */
void casimir_free(casimir_t *self)
{
    //gsl_integration_workspace_free(self->int_workspace);
}


/*
 * For small x<<1 a_l will scale as
 * 
 * a_l(x) ~ a0*(x/2)^(2l+1)
 *
 * This method returns the prefactor a0
 */
double casimir_a0(int l)
{
    return M_PI*pow(-1, l)*( 2*gsl_sf_gamma(l+1.5)-l*gsl_sf_gamma(l+0.5) )/( l*pow(gsl_sf_gamma(l+0.5),2)*gsl_sf_gamma(l+1.5) );
}

/*
 * For small x<<1 b_l will scale as
 * 
 * b_l(x) ~ b0*(x/2)^(2l+1)
 *
 * This method returns the prefactor b0
 */
double casimir_b0(int l)
{
    return M_PI*pow(-1, l+1)/( gsl_sf_gamma(l+0.5)*gsl_sf_gamma(l+1.5) );
}

/*
 * Returns the coefficient a_l for reflection on the sphere
 *
 * Restrictions: l integer, l>=1, xi>0
 */
double inline casimir_a(casimir_t *self, int l, double arg)
{
    return M_PI/2*pow(-1,l+1)*(l*gsl_sf_bessel_Inu(l+0.5,arg)-arg*gsl_sf_bessel_Inu(l-0.5,arg))/(l*gsl_sf_bessel_Knu(l+0.5,arg)+arg*gsl_sf_bessel_Knu(l-0.5,arg));
}

/*        
 * Returns the coefficient b_l for reflection on the sphere
 *
 * Restrictions: l integer, l>=1, xi>0
 */        
double casimir_b(casimir_t *self, int  l, double arg)
{
    return M_PI/2*pow(-1, l+1)*gsl_sf_bessel_Inu(l+0.5,arg)/gsl_sf_bessel_Knu(l+0.5,arg);
}

/*
 * Returns the TE reflection coefficient for the plate
 */
inline double casimir_rTE(casimir_t *self, double x, double arg)
{
    return -1;
}

/*
 * Returns the TM reflection coefficient for the plate
 */
double casimir_rTM(casimir_t *self, double x, double xi)
{
    return +1;
}

/*
 * Initialize the mie cache.
 * This function must be called before any call to casimir_mie_cache_alloc
 */
void casimir_mie_cache_init(casimir_mie_cache_t *cache, double arg)
{
    cache->al = cache->bl = NULL;
    cache->lmax = 0;
    cache->arg  = arg;
}

/*
 * Allocate memory for the Mie-coefficients a_l and b_l
 */
int casimir_mie_cache_alloc(casimir_t *self, casimir_mie_cache_t *cache, int lmax)
{
    int l1;
    double arg = cache->arg;

    if(arg == 0)
    {
        cache->al = cache->bl = NULL;
        return 0;
    }

    cache->al = (double *)realloc(cache->al, (lmax+1)*sizeof(double));
    cache->bl = (double *)realloc(cache->bl, (lmax+1)*sizeof(double));

    if(cache->al == NULL || cache->bl == NULL)
    {
        fprintf(stderr, "# Out of memory.\n");
        exit(1);
    }

    cache->al[0] = cache->bl[0] = 0;
    for(l1 = MAX(1,cache->lmax); l1 <= lmax; l1++)
    {
        cache->al[l1] = casimir_a(self,l1,arg);
        cache->bl[l1] = casimir_b(self,l1,arg);
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
        free(cache->al);
    if(cache->bl != NULL)
        free(cache->bl);

    cache->al = cache->bl = NULL;
}

/*
 * Calculate free energy. Sum from nmin to nmax.

 * Restrictions: nmax integer, nmax >= 0
 */
double casimir_F(casimir_t *self, int *nmax)
{
    int n = 0;
    double sum_n0 = 0, sum = 0;
    double RbyScriptL = self->RbyScriptL;

    while(1)
    {
        casimir_mie_cache_t cache;
        double sum_n = 0;
        int m;

        casimir_mie_cache_init(&cache, n*self->T*RbyScriptL);
        casimir_mie_cache_alloc(self, &cache, self->lmax);

        for(m = 0; m <= self->lmax; m++)
        {
            double value = casimir_logdetD(self,n,m,&cache);

            if(m == 0)
                value /= 2;

            sum_n += value;
        }

        if(n == 0)
        {
            sum_n0 = sum_n;
            sum += sum_n/2;
        }
        else
            sum += sum_n;

        //fprintf(stderr, "# n=%d, %.15g\n", n, sum_n);
        casimir_mie_cache_free(&cache);

        if(sum_n/sum_n0 < self->eps_n)
        {
            if(nmax != NULL)
                *nmax = n;
            return self->T/M_PI*sum;
        }

        n++;
    }
}

/*
 * Calculate logarithm of determinant of D=1-M for given m,ξ
 *
 * Restrictions: m integer, m>=0, ξ>= 0
 */
double casimir_logdetD(casimir_t *self, int n, int m, casimir_mie_cache_t *cache)
{
    int min,max,dim,l1,l2;
    double logdet_EE = 0;
    double logdet_MM = 0;
    double logdet = 0;

    min = MAX(m,1);
    max = self->lmax;
    
    dim = (max-min+1);
    
    if(n == 0)
    {
        gsl_matrix *EE, *MM;
        double RbyScriptL = self->RbyScriptL;

        EE = gsl_matrix_alloc(dim, dim);
        MM = gsl_matrix_alloc(dim, dim);

        for(l1 = min; l1 <= max; l1++)
            for(l2 = min; l2 <= max; l2++)
            {
                double XiRL = casimir_Xi(l1,l2,m)*pow(RbyScriptL, 2*l1+1);

                gsl_matrix_set(EE, l1-min, l2-min, +casimir_a0(l1)*XiRL);
                gsl_matrix_set(MM, l1-min, l2-min, -casimir_b0(l1)*XiRL);
            }
    
        casimir_logdet1m(EE, &logdet_EE, n, m, "M_EE");
        casimir_logdet1m(MM, &logdet_MM, n, m, "M_MM");

        gsl_matrix_free(EE);
        gsl_matrix_free(MM);

        return logdet_EE+logdet_MM;
    }
    else
    {
        //double nTRbyScriptL = n*self->T*self->RbyScriptL;
        gsl_matrix *M = gsl_matrix_alloc(2*dim, 2*dim);
    
        /* M_EE, -M_EM
           M_ME,  M_MM */
        for(l1 = min; l1 <= max; l1++)
        {
            for(l2 = min; l2 <= l1; l2++)
            {
                casimir_integrals_t cint;
                double scale = 1;
                double al1 = cache->al[l1];
                double bl1 = cache->bl[l1];
                double al2 = cache->al[l2];
                double bl2 = cache->bl[l2];

                /*
                if(nTRbyScriptL < 1)
                {
                    double scale_l1, scale_l2;
                    scale    = pow(nTRbyScriptL, +l2+l1);
                    scale_l1 = pow(nTRbyScriptL, 2*l1);
                    scale_l2 = pow(nTRbyScriptL, 2*l2);

                    al1 /= scale_l1;
                    bl1 /= scale_l1;
                    al2 /= scale_l2;
                    bl2 /= scale_l2;

                    if(fabs(al1) < 1e-250)
                        al1 = casimir_a0(l1)*nTRbyScriptL/pow(2, 2*l1+1);
                    if(fabs(bl1) < 1e-250)
                        bl1 = casimir_b0(l1)*nTRbyScriptL/pow(2, 2*l1+1);
                    if(fabs(al2) < 1e-250)
                        al2 = casimir_a0(l2)*nTRbyScriptL/pow(2, 2*l2+1);
                    if(fabs(bl2) < 1e-250)
                        bl2 = casimir_b0(l2)*nTRbyScriptL/pow(2, 2*l2+1);
                }
                */

                casimir_integrate(&cint, l1, l2, m, 2*n*self->T, scale);

                gsl_matrix_set(M,     l1-min,     l2-min, -2*al1*(cint.B-cint.A)); /* M_EE */
                gsl_matrix_set(M,     l2-min,     l1-min, -2*pow(-1, l1+l2)*al2*(cint.B-cint.A)); /* M_EE */

                gsl_matrix_set(M, dim+l1-min, dim+l2-min, -2*bl1*(cint.A-cint.B)); /* M_MM */
                gsl_matrix_set(M, dim+l2-min, dim+l1-min, -2*pow(-1, l1+l2)*bl2*(cint.A-cint.B)); /* M_MM */

                gsl_matrix_set(M, dim+l1-min,     l2-min, -2*al1*(cint.D-cint.C)); /* M_EM */
                gsl_matrix_set(M, dim+l2-min,     l1-min, -2*pow(-1, l1+l2+1)*al2*(cint.D-cint.C)); /* M_EM */

                gsl_matrix_set(M,     l1-min, dim+l2-min, -2*bl1*(cint.C-cint.D)); /* - M_ME */
                gsl_matrix_set(M,     l2-min, dim+l1-min, -2*pow(-1, l1+l2+1)*bl2*(cint.C+cint.D)); /* - M_ME */
            }
        }

        /*
        fprintf(stderr, "%d\n", (int)M->size1);
        gsl_matrix_fprintf(stderr, M, "%g");
        */

        if(m == 0)
        {
            size_t i,j;
            gsl_matrix *EE = gsl_matrix_alloc(dim, dim);
            gsl_matrix *MM = gsl_matrix_alloc(dim, dim);

            for(i = 0; i < dim; i++)
                for(j = 0; j < dim; j++)
                {
                    gsl_matrix_set(EE, i,j, gsl_matrix_get(M, i,j));
                    gsl_matrix_set(MM, i,j, gsl_matrix_get(M, dim+i,dim+j));
                }

            casimir_logdet1m(EE, &logdet_EE, n, m, "EE");
            casimir_logdet1m(MM, &logdet_MM, n, m, "MM");

            logdet = logdet_EE+logdet_MM;

            gsl_matrix_free(EE);
            gsl_matrix_free(MM);
        }
        else
            casimir_logdet1m(M, &logdet, n, m, "M");

        gsl_matrix_free(M);
        return logdet;
    }
}
