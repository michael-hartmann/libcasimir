#include <math.h>
#include <assert.h>
#include <quadmath.h>

#include "givens.h"
#include "sfunc.h"
#include "integration.h"
#include "libcasimir.h"

#define FACTOR_LMAX 5
#define HBARC 3.161526510740123e-26
#define KB   1.3806488e-23

#define EPS_PRECISION 1e-16

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
    return (log(2.*l1+1)+log(2*l2+1)-log(4)-log(l1)-log(l1+1)-log(l2)-log(l2+1)+lnfac(l1-m)+lnfac(l2-m)-lnfac(l1+m)-lnfac(l2+m))/2.0;
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
           +lnfac(2*l1)+lnfac(2*l2)+lnfac(l1+l2)-log(4)*(2*l1+l2+1)-lnfac(l1-1)-lnfac(l2-1);
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

    self->T          = T;
    self->RbyScriptL = RbyScriptL;
    self->precision  = EPS_PRECISION;
    self->verbose    = 0;

    return 0;
}

/*
 * Set maximum value for l. lmax must be positive.
 */
void casimir_set_lmax(casimir_t *self, int lmax)
{
    if(lmax > 0)
        self->lmax = lmax;
}

double casimir_get_lmax(casimir_t *self)
{
    return self->lmax;
}

void casimir_set_verbose(casimir_t *self, int verbose)
{
    self->verbose = verbose ? 1 : 0;
}

void casimir_set_precision(casimir_t *self, double precision)
{
    if(precision > 0)
        self->precision = precision;
}

/*
 * Free memory for casimir object
 */
void casimir_free(casimir_t *self)
{
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

    /* we could do both calculation together. but it doesn't cost much time -
     * so why bother? 
     */
    bessel_lnInuKnu(l-1, arg, &lnIlm, &lnKlm);
    bessel_lnInuKnu(l,   arg, &lnIlp, &lnKlp);

    prefactor = M_LOGPI-log(2)+lnIlp-lnKlp;
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

    cache->al      = (double *)realloc(cache->al,      (lmax+1)*sizeof(double));
    cache->bl      = (double *)realloc(cache->bl,      (lmax+1)*sizeof(double));
    cache->al_sign =    (int *)realloc(cache->al_sign, (lmax+1)*sizeof(int));
    cache->bl_sign =    (int *)realloc(cache->bl_sign, (lmax+1)*sizeof(int));

    if(cache->al == NULL || cache->bl == NULL || cache->al_sign == NULL || cache->bl_sign == NULL)
    {
        fprintf(stderr, "# Out of memory.\n");
        exit(1);
    }

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
        free(cache->al);
    if(cache->bl != NULL)
        free(cache->bl);
    if(cache->al_sign != NULL)
        free(cache->al_sign);
    if(cache->bl_sign != NULL)
        free(cache->bl_sign);

    cache->al = cache->bl = NULL;
}

/*
 * Calculate free energy.

 * Restrictions: nmax integer, nmax >= 0
 */
double casimir_F(casimir_t *self, int *nmax)
{
    int n = 0;
    double F = 0;
    double precision = self->precision;
    double TRbyScriptL = self->T*self->RbyScriptL;

    /* So, here we sum up all m and n that contribute to F.
     * So, what do we do here?
     *
     * We want to evaluate
     *      \sum_{n=0}^\infty \sum{m=0}^{l_max} log(det(1-M)),
     * where the terms for n=0 and m=0 are weighted with a factor 1/2.
     */
    while(1)
    {
        casimir_mie_cache_t cache;
        int m;

        casimir_mie_cache_init(&cache, n*TRbyScriptL);
        casimir_mie_cache_alloc(self, &cache, self->lmax);

        double sum_n = 0;
        for(m = 0; m <= self->lmax; m++)
        {
            double value = casimir_logdetD(self,n,m,&cache);
            if(self->verbose)
                fprintf(stderr, "# n=%d, m=%d, value=%.15g\n", n, m, value);

            /* If F is !=0 and value/F < 1e-16, then F+value = F. The addition
             * has no effect.
             * As for larger m value will be even smaller, we can skip the
             * summation here. 
             */
            if(F != 0 && fabs(value/F) < precision)
            {
                F += value;
                break;
            }
            if(m == 0)
                value /= 2;
            if(n == 0)
                value /= 2;

            F += value;
            sum_n += value;
        }
        fprintf(stderr, "# n=%d, value=%.15g\n", n, sum_n);

        casimir_mie_cache_free(&cache);

        /* If m == 0, all other terms will be even smaller and we can skip the
         * summation.
         */
        if(m == 0)
        {
            /* get out of here */
            if(nmax != NULL)
                *nmax = n;
            return self->T/M_PI*F;
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
        double RbyScriptL = self->RbyScriptL;
        matrix_t *EE, *MM;
        EE = matrix_alloc(dim);
        MM = matrix_alloc(dim);

        for(l1 = min; l1 <= max; l1++)
            for(l2 = min; l2 <= max; l2++)
            {
                int sign_a0, sign_b0, sign_xi;
                __float128 elem_EE, elem_MM;
                double ln_a0, ln_b0;
                double lnXiRL = casimir_lnXi(l1,l2,m,&sign_xi)+(2*l1+1)*log(RbyScriptL);
                casimir_lna0_lnb0(l1, &ln_a0, &sign_a0, &ln_b0, &sign_b0);

                elem_EE = (l1 == l2 ? 1 : 0)-(sign_xi*sign_a0)*expq(ln_a0+lnXiRL);
                elem_MM = (l1 == l2 ? 1 : 0)+(sign_xi*sign_b0)*expq(ln_b0+lnXiRL);

                assert(!isnanq(elem_EE));
                assert(!isinfq(elem_EE));
                assert(!isnanq(elem_MM));
                assert(!isinfq(elem_MM));

                matrix_set(EE, l1-min, l2-min, elem_EE);
                matrix_set(MM, l1-min, l2-min, elem_MM);
            }
    
        logdet_EE = matrix_logdet(EE);
        logdet_MM = matrix_logdet(MM);

        matrix_free(EE);
        matrix_free(MM);

        return logdet_EE+logdet_MM;
    }
    else
    {
        double nTRbyScriptL = n*self->T*self->RbyScriptL;
        matrix_t *M = matrix_alloc(2*dim);
    
        /* M_EE, -M_EM
           M_ME,  M_MM */
        for(l1 = min; l1 <= max; l1++)
        {
            for(l2 = min; l2 <= l1; l2++)
            {
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

                casimir_integrate(&cint, l1, l2, m, 2*n*self->T);

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

                /* M_EE */
                {
                    __float128 M_EE1 = -2*al1_sign               *( cint.signB*expq(lnal1+cint.logB)-cint.signA*expq(lnal1+cint.logA) );
                    __float128 M_EE2 = -2*al2_sign*pow(-1, l1+l2)*( cint.signB*expq(lnal2+cint.logB)-cint.signA*expq(lnal2+cint.logA) );

                    assert(!isnanq(M_EE1)); assert(!isinfq(M_EE1));
                    assert(!isnanq(M_EE2)); assert(!isinfq(M_EE2));

                    matrix_set(M, l1-min, l2-min, (l1 == l2 ? 1 : 0)-M_EE1); /* M_EE */
                    matrix_set(M, l2-min, l1-min, (l1 == l2 ? 1 : 0)-M_EE2); /* M_EE */
                }

                /* M_MM */
                {
                    __float128 M_MM1 = -2*bl1_sign               *( cint.signA*expq(lnbl1+cint.logA)-cint.signB*expq(lnbl1+cint.logB) );
                    __float128 M_MM2 = -2*bl2_sign*pow(-1, l1+l2)*( cint.signA*expq(lnbl2+cint.logA)-cint.signB*expq(lnbl2+cint.logB) );

                    assert(!isnanq(M_MM1)); assert(!isinfq(M_MM1));
                    assert(!isnanq(M_MM2)); assert(!isinfq(M_MM2));

                    matrix_set(M, dim+l1-min, dim+l2-min, (l1 == l2 ? 1 : 0)-M_MM1); /* M_MM */
                    matrix_set(M, dim+l2-min, dim+l1-min, (l1 == l2 ? 1 : 0)-M_MM2); /* M_MM */
                }

                /* M_EM */
                {
                    __float128 M_EM1 = -2*al1_sign                *( cint.signD*expq(lnal1+cint.logD)-cint.signC*expq(lnal1+cint.logC) );
                    __float128 M_EM2 = -2*al2_sign*pow(-1,l1+l2+1)*( cint.signC*expq(lnal2+cint.logC)-cint.signD*expq(lnal2+cint.logD) );

                    assert(!isnanq(M_EM1)); assert(!isinfq(M_EM1));
                    assert(!isnanq(M_EM2)); assert(!isinfq(M_EM2));

                    matrix_set(M, dim+l1-min, l2-min, -M_EM1); /* M_EM */
                    matrix_set(M, dim+l2-min, l1-min, -M_EM2); /* M_EM */
                }


                /* M_ME */
                {
                    __float128 M_ME1 = -2*bl1_sign                *( cint.signC*expq(lnbl1+cint.logC)-cint.signD*expq(lnbl1+cint.logD) );
                    __float128 M_ME2 = -2*bl2_sign*pow(-1,l1+l2+1)*( cint.signD*expq(lnbl2+cint.logD)-cint.signC*expq(lnbl2+cint.logC) );

                    assert(!isnanq(M_ME1)); assert(!isinfq(M_ME1));
                    assert(!isnanq(M_ME2)); assert(!isinfq(M_ME2));

                    matrix_set(M, l1-min, dim+l2-min, -M_ME1); /* - M_ME */
                    matrix_set(M, l2-min, dim+l1-min, -M_ME2); /* - M_ME */
                }
            }
        }

        /*
        fprintf(stderr, "%d\n", (int)M->size1);
        gsl_matrix_fprintf(stderr, M, "%g");
        */

        if(m == 0)
        {
            size_t i,j;
            matrix_t *EE = matrix_alloc(dim);
            matrix_t *MM = matrix_alloc(dim);

            for(i = 0; i < dim; i++)
                for(j = 0; j < dim; j++)
                {
                    matrix_set(EE, i,j, matrix_get(M, i,j));
                    matrix_set(MM, i,j, matrix_get(M, dim+i,dim+j));
                }

            logdet = matrix_logdet(EE)+matrix_logdet(MM);

            matrix_free(EE);
            matrix_free(MM);
        }
        else
            logdet = matrix_logdet(M);


        matrix_free(M);
        assert(!isinf(logdet));
        return logdet;
    }
}
