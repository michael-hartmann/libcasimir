#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "matrix.h"
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

        if(self->verbose)
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
        double lnRbyScriptL = log(self->RbyScriptL);

        #ifdef MATRIX_QUAD
            matrix_quad_t *EE = matrix_quad_alloc(dim);
            matrix_quad_t *MM = matrix_quad_alloc(dim);
        #else
            matrix_char_t *EE_signs = matrix_char_alloc(dim);
            matrix_char_t *MM_signs = matrix_char_alloc(dim);
            matrix_t *EE       = matrix_alloc(dim);
            matrix_t *MM       = matrix_alloc(dim);
        #endif

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

                #ifdef MATRIX_QUAD
                    matrix_set(EE, i,j, (l1 == l2 ? 1 : 0) - sign_xi*sign_a0*expq(lna0+lnXiRL));
                    matrix_set(MM, i,j, (l1 == l2 ? 1 : 0) - sign_xi*sign_a0*expq(lnb0+lnXiRL));
                #else
                    double D_EE, D_MM;
                    int sign_EE, sign_MM;

                    /* -M_EE */
                    D_EE    = lna0+lnXiRL;
                    sign_EE = -(sign_xi*sign_a0);

                    /* -M_MM */
                    D_MM    = lnb0+lnXiRL;
                    sign_MM = +(sign_xi*sign_b0);

                    /* if l1 == l2: D_EE = 1-M_EE and D_MM = 1-M_MM */
                    if(l1 == l2)
                    {
                        D_EE = logadd_s(0, +1, D_EE, sign_EE, &sign_EE);
                        D_MM = logadd_s(0, +1, D_MM, sign_MM, &sign_MM);
                    }

                    matrix_set(EE_signs, i,j, sign_EE);
                    matrix_set(EE, i, j, D_EE);

                    matrix_set(MM_signs, i,j, sign_MM);
                    matrix_set(MM, i, j, D_MM);
                #endif
            }
    
        #ifdef MATRIX_QUAD
            /* balance the matrix */
            matrix_quad_balance(EE);
            matrix_quad_balance(MM);

            logdet_EE = matrix_quad_logdet(EE);
            logdet_MM = matrix_quad_logdet(MM);

            /* free space for matrices */
            matrix_quad_free(EE);
            matrix_quad_free(MM);
        #else
            /* balance the matrix */
            matrix_log_balance(EE);
            matrix_log_balance(MM);

            /* exp the matrix */
            matrix_exp(EE);
            matrix_exp(MM);

            /* multiply the correct signs */
            for(l1 = 0; l1 < dim; l1++)
                for(l2 = 0; l2 < dim; l2++)
                {
                    matrix_set(EE, l1,l2, matrix_get(EE_signs, l1,l2)*matrix_get(EE,l1,l2));
                    matrix_set(MM, l1,l2, matrix_get(MM_signs, l1,l2)*matrix_get(MM,l1,l2));
                }

            /* free space for sign matrices */
            matrix_char_free(EE_signs);
            matrix_char_free(MM_signs);

            /* now EE and MM contain the "normal" matrix elements, so we can simply
             * calculate the determinant of EE and MM 
             */
            logdet_EE = matrix_logdet(EE);
            logdet_MM = matrix_logdet(MM);

            /* free space for matrices */
            matrix_free(EE);
            matrix_free(MM);
        #endif

        return logdet_EE+logdet_MM;
    }
    else
    {
        double nTRbyScriptL = n*self->T*self->RbyScriptL;
        #ifdef MATRIX_QUAD
            matrix_quad_t *M = matrix_quad_alloc(2*dim);
        #else
            matrix_t *M = matrix_alloc(2*dim);
            matrix_char_t *M_signs = matrix_char_alloc(2*dim);
        #endif
    
        /* M_EE, -M_EM
           M_ME,  M_MM */
        for(l1 = min; l1 <= max; l1++)
        {
            for(l2 = min; l2 <= l1; l2++)
            {
                const int i = l1-min, j = l2-min;
                #ifndef MATRIX_QUAD
                    int sign;
                #endif
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
                #ifdef MATRIX_QUAD
                    /* EE */
                    matrix_set(M, i,j, (l1 == l2 ? 1 : 0) + 2*al1_sign*( cint.signB*expq(lnal1+cint.logB) - cint.signA*expq(lnal1+cint.logA) ));
                    matrix_set(M, j,i, (l1 == l2 ? 1 : 0) + 2*pow(-1, l1+l2)*al2_sign*( cint.signB*expq(lnal2+cint.logB) - cint.signA*expq(lnal2+cint.logA) ));

                    /* MM */
                    matrix_set(M, i+dim,j+dim, (l1 == l2 ? 1 : 0) + 2*bl1_sign*( cint.signA*expq(lnbl1+cint.logA) - cint.signB*expq(lnbl1+cint.logB) ));
                    matrix_set(M, j+dim,i+dim, (l1 == l2 ? 1 : 0) + 2*pow(-1, l1+l2)*bl2_sign*( cint.signA*expq(lnbl2+cint.logA) - cint.signB*expq(lnbl2+cint.logB) ));
                #else
                    double D_EE1 = M_LN2 + logadd_s(lnal1+cint.logB, cint.signB, lnal1+cint.logA, -cint.signA, &sign);
                    int sign_EE1 = sign*al1_sign;
                    double D_EE2 = M_LN2 + logadd_s(lnal2+cint.logB, cint.signB, lnal2+cint.logA, -cint.signA, &sign);
                    int sign_EE2 = sign*al2_sign*pow(-1, l1+l2);

                    double D_MM1 = M_LN2 + logadd_s(lnbl1+cint.logA, cint.signA, lnbl1+cint.logB, -cint.signB, &sign);
                    int sign_MM1 = sign*bl1_sign;
                    double D_MM2 = M_LN2 + logadd_s(lnbl2+cint.logA, cint.signA, lnbl2+cint.logB, -cint.signB, &sign);
                    int sign_MM2 = sign*bl2_sign*pow(-1, l1+l2);

                    if(l1 == l2)
                    {
                        D_EE1 = logadd_s(0, +1, D_EE1, sign_EE1, &sign_EE1);
                        D_EE2 = logadd_s(0, +1, D_EE2, sign_EE2, &sign_EE2);

                        D_MM1 = logadd_s(0, +1, D_MM1, sign_MM1, &sign_MM1);
                        D_MM2 = logadd_s(0, +1, D_MM2, sign_MM2, &sign_MM2);
                    }

                    assert(!isnan(D_EE1)); assert(!isinf(D_EE1));
                    assert(!isnan(D_EE2)); assert(!isinf(D_EE2));
                    assert(!isnan(D_MM1)); assert(!isinf(D_MM1));
                    assert(!isnan(D_MM2)); assert(!isinf(D_MM2));

                    matrix_set(M_signs, i,j, sign_EE1);
                    matrix_set(M, i,j, D_EE1); /* M_EE */

                    matrix_set(M_signs, j,i, sign_EE2);
                    matrix_set(M, j,i, D_EE2);

                    matrix_set(M_signs, i+dim,j+dim, sign_MM1);
                    matrix_set(M, i+dim,j+dim, D_MM1); /* M_EE */

                    matrix_set(M_signs, j+dim,i+dim, sign_MM2);
                    matrix_set(M, j+dim,i+dim, D_MM2);
                #endif


                if(m != 0)
                {
                    /* M_EM */
                    #ifdef MATRIX_QUAD
                        matrix_set(M, dim+i,j, 2*al1_sign*( cint.signD*expq(lnal1+cint.logD) - cint.signC*expq(lnal1+cint.logC) ));
                        matrix_set(M, dim+j,i, 2*al2_sign*pow(-1, l1+l2+1)*( cint.signC*expq(lnal2+cint.logC) - cint.signD*expq(lnal2+cint.logD) ));
                    #else
                        double D_EM1 = M_LN2 + logadd_s(lnal1+cint.logD, cint.signD, lnal1+cint.logC, -cint.signC, &sign);
                        int sign_EM1 = sign*al1_sign;
                        double D_EM2 = M_LN2 + logadd_s(lnal2+cint.logC, cint.signC, lnal2+cint.logD, -cint.signD, &sign);
                        int sign_EM2 = sign*al2_sign*pow(-1, l1+l2+1);

                        assert(!isnan(D_EM1)); assert(!isinf(D_EM1));
                        assert(!isnan(D_EM2)); assert(!isinf(D_EM2));

                        matrix_set(M_signs, dim+i,j, sign_EM1);
                        matrix_set(M, dim+i, j, D_EM1); /* M_EM */

                        matrix_set(M_signs, dim+j,i, sign_EM2);
                        matrix_set(M, dim+j, i, D_EM2); /* M_EM */
                    #endif


                    /* M_ME */
                    #ifdef MATRIX_QUAD
                        matrix_set(M, i,dim+j, 2*bl1_sign*( cint.signC*expq(lnbl1+cint.logC) - cint.signD*expq(lnbl1+cint.logD) ));
                        matrix_set(M, j,dim+i, 2*bl2_sign*pow(-1, l1+l2+1)*( cint.signD*expq(lnbl2+cint.logD) - cint.signC*expq(lnbl2+cint.logC) ));
                    #else
                        double D_ME1 = M_LN2 + logadd_s(lnbl1+cint.logC, cint.signC, lnbl1+cint.logD, -cint.signD, &sign);
                        int sign_ME1 = sign*bl1_sign;
                        double D_ME2 = M_LN2 + logadd_s(lnbl2+cint.logD, cint.signD, lnbl2+cint.logC, -cint.signC, &sign);
                        int sign_ME2 = sign*bl2_sign*pow(-1, l1+l2+1);

                        assert(!isnan(D_ME1)); assert(!isinf(D_ME1));
                        assert(!isnan(D_ME2)); assert(!isinf(D_ME2));

                        matrix_set(M_signs, i,dim+j, sign_ME1);
                        matrix_set(M, i, dim+j, D_ME1); /* - M_ME */

                        matrix_set(M_signs, j,dim+i, sign_ME2);
                        matrix_set(M, j, dim+i, D_ME2); /* - M_ME */
                    #endif
                }
            }
        }

        if(m == 0)
        {
            size_t i,j;
            #ifdef MATRIX_QUAD
                matrix_quad_t *EE = matrix_quad_alloc(dim);
                matrix_quad_t *MM = matrix_quad_alloc(dim);
            #else
                matrix_t *EE = matrix_alloc(dim);
                matrix_t *MM = matrix_alloc(dim);
            #endif

            for(i = 0; i < dim; i++)
                for(j = 0; j < dim; j++)
                {
                    matrix_set(EE, i,j, matrix_get(M, i,j));
                    matrix_set(MM, i,j, matrix_get(M, dim+i,dim+j));
                }

            #ifdef MATRIX_QUAD
                matrix_quad_balance(MM);
                matrix_quad_balance(EE);

                logdet = matrix_quad_logdet(EE)+matrix_quad_logdet(MM);

                matrix_quad_free(EE);
                matrix_quad_free(MM);
            #else
                matrix_log_balance(EE);
                matrix_exp(EE);

                matrix_log_balance(MM);
                matrix_exp(MM);

                for(l1 = 0; l1 < dim; l1++)
                    for(l2 = 0; l2 < dim; l2++)
                    {
                        matrix_set(EE, l1,l2, matrix_get(M_signs, l1,l2)        *matrix_get(EE,l1,l2));
                        matrix_set(MM, l1,l2, matrix_get(M_signs, l1+dim,l2+dim)*matrix_get(MM,l1,l2));
                    }

                /* free space for sign matrices */
                matrix_char_free(M_signs);

                /* calculate log det */
                logdet = matrix_logdet(EE)+matrix_logdet(MM);

                matrix_free(EE);
                matrix_free(MM);
            #endif
        }
        else
        {
            #ifdef MATRIX_QUAD
                matrix_quad_balance(M);
                logdet = matrix_quad_logdet(M);
            #else
                matrix_log_balance(M);
                matrix_exp(M);

                for(l1 = 0; l1 < 2*dim; l1++)
                    for(l2 = 0; l2 < 2*dim; l2++)
                        matrix_set(M, l1,l2, matrix_get(M_signs, l1,l2)*matrix_get(M,l1,l2));

                /* free space for sign matrices */
                matrix_char_free(M_signs);

                /* calculate log det */
                logdet = matrix_logdet(M);
            #endif
        }

        #ifdef MATRIX_QUAD
            matrix_quad_free(M);
        #else
            matrix_free(M);
        #endif

        assert(!isinf(logdet));
        return logdet;
    }
}
