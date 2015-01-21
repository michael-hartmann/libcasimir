#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "utils.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration.h"

/* weights for Gauss-Legendre integration order 25 */
double ln_wk[] = {
    -1.779355901402488,
    -1.233558297369509,
    -1.321682956887102,
    -1.795752712002829,
    -2.592588996560207,
    -3.690303775607052,
    -5.082794644999876,
    -6.772383492754023,
    -8.767311223343112,
    -11.08104949024446,
    -13.7325615854023,
    -16.74730682131207,
    -20.15909974362777,
    -24.01323430703988,
    -28.37177801482692,
    -33.32302001737665,
    -38.99983382185022,
    -45.62022978993102,
    -53.596896664818,
    -63.96770185369194
};


/* nodes for Gauss-Legendre integration order 25 */
double xk[] = {
    0.705398896919887533667e-1,
    0.372126818001611443794,
    0.916582102483273564668,
    0.170730653102834388069e1,
    0.274919925530943212965e1,
    0.404892531385088692237e1,
    0.561517497086161651410e1,
    0.745901745367106330977e1,
    0.959439286958109677247e1,
    0.120388025469643163096e2,
    0.148142934426307399785e2,
    0.179488955205193760174e2,
    0.214787882402850109757e2,
    0.254517027931869055035e2,
    0.299325546317006120067e2,
    0.350134342404790000063e2,
    0.408330570567285710620e2,
    0.476199940473465021399e2,
    0.558107957500638988908e2,
    0.665244165256157538186e2
};

void integrands_drude(double x, integrands_drude_t *integrands, casimir_t *self, double nT, int l1, int l2, int m)
{
    plm_combination_t comb;
    const double tau = 2*nT;
    const double k = sqrt(pow_2(x)/4 + nT*x);
    const double log_factor = log(pow_2(x)+2*tau*x);
    double r_TE, r_TM, lnr_TE, lnr_TM;
    double A,B,C,D;

    casimir_rp(self, nT, k, &r_TE, &r_TM);
    lnr_TE  = log(-r_TE);
    lnr_TM = log(r_TM);

    plm_PlmPlm(l1, l2, m, 1+x/tau, &comb);

    A = comb.lnPl1mPl2m - log_factor;
    integrands->lnA_TE = lnr_TE + A;
    integrands->lnA_TM = lnr_TM + A;
    integrands->sign_A = comb.sign_Pl1mPl2m;

    B = comb.lndPl1mdPl2m + log_factor;
    integrands->lnB_TE = lnr_TE + B;
    integrands->lnB_TM = lnr_TM + B;
    integrands->sign_B = comb.sign_dPl1mdPl2m;

    C = comb.lnPl1mdPl2m;
    integrands->lnC_TE = lnr_TE + C;
    integrands->lnC_TM = lnr_TM + C;
    integrands->sign_C = comb.sign_Pl1mdPl2m;

    D = comb.lndPl1mPl2m;
    integrands->lnD_TE = lnr_TE + D;
    integrands->lnD_TM = lnr_TM + D;
    integrands->sign_D = comb.sign_dPl1mPl2m;
}


/** @brief Calculate integrals A,B,C,D including prefactor Lambda vor Drude metals
 *
 * This function calculates
 *    Lambda(l1,l2,m)*A_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*B_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*C_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*D_(l1,l2)^(m)
 * for Drude metals.
 *
 * @param [in]  self Casimir object
 * @param [out] cint logarithms of values and signs of integrals
 * @param [in]  l1   \f$\ell_1\f$
 * @param [in]  l2   \f$\ell_2\f$
 * @param [in]  m    \f$m\f$
 * @param [in]  nT   \f$nT\f$
 */
void casimir_integrate_drude(casimir_t *self, casimir_integrals_t *cint, int l1, int l2, int m, double nT)
{
    int i;
    const int N = sizeof(ln_wk)/sizeof(double);
    integrands_drude_t integrand;
    const double tau = 2*nT;
    const double ln_tau = log(2*nT);
    const double ln_Lambda = casimir_lnLambda(l1, l2, m, NULL); /* sign: -1 */
    double prefactor;
    double *ln_ABCD, *lnA_TE, *lnA_TM, *lnB_TE, *lnB_TM, *lnC_TE, *lnC_TM, *lnD_TE, *lnD_TM;
    int *signs_ABCD, *signs_A, *signs_B, *signs_C, *signs_D;

    /* allocate space for signs_A, signs_B, signs_C, signs_D */
    signs_ABCD = xmalloc(4*N*sizeof(int));
    signs_A = signs_ABCD;
    signs_B = signs_A+1*N;
    signs_C = signs_A+2*N;
    signs_D = signs_A+3*N;

    /* allocate space for lnA_TE, lnA_TM, lnB_TE, lnB_TM, lnC_TE, lnC_TM,
     * lnD_TE, lnD_TM */
    ln_ABCD = xmalloc(4*2*N*sizeof(integrands_drude_t));
    lnA_TE  = ln_ABCD;
    lnA_TM  = ln_ABCD + 1*N;
    lnB_TE  = ln_ABCD + 2*N;
    lnB_TM  = ln_ABCD + 3*N;
    lnC_TE  = ln_ABCD + 4*N;
    lnC_TM  = ln_ABCD + 5*N;
    lnD_TE  = ln_ABCD + 6*N;
    lnD_TM  = ln_ABCD + 7*N;

    for(i = 0; i < N; i++)
    {
        integrands_drude(xk[i], &integrand, self, nT, l1, l2, m);

        lnA_TE[i]  = ln_wk[i] + integrand.lnA_TE;
        lnA_TM[i]  = ln_wk[i] + integrand.lnA_TM;
        signs_A[i] = integrand.sign_A;

        lnB_TE[i]  = ln_wk[i] + integrand.lnB_TE;
        lnB_TM[i]  = ln_wk[i] + integrand.lnB_TM;
        signs_B[i] = integrand.sign_B;

        lnC_TE[i]  = ln_wk[i] + integrand.lnC_TE;
        lnC_TM[i]  = ln_wk[i] + integrand.lnC_TM;
        signs_C[i] = integrand.sign_C;

        lnD_TE[i]  = ln_wk[i] + integrand.lnD_TE;
        lnD_TM[i]  = ln_wk[i] + integrand.lnD_TM;
        signs_D[i] = integrand.sign_D;
    }


    /* B */
    prefactor = -tau-3*ln_tau; /* exp(-tau)/tau³ */
    cint->lnB_TE = prefactor + logadd_ms(lnB_TE, signs_B, N, &cint->signB_TE);
    cint->lnB_TM = prefactor + logadd_ms(lnB_TM, signs_B, N, &cint->signB_TM);

    cint->signB_TM = -pow(-1, l2+m+1) * cint->signB_TM;
    cint->signB_TE = +pow(-1, l2+m+1) * cint->signB_TE;


    if(m > 0)
    {
        const double log_m = log(m);

        /* A */
        prefactor = 2*log_m+ln_tau-tau+ln_Lambda; /* m²*tau*exp(-tau) */
        cint->lnA_TE = prefactor + logadd_ms(lnA_TE, signs_A, N, &cint->signA_TE);
        cint->lnA_TM = prefactor + logadd_ms(lnA_TM, signs_A, N, &cint->signA_TM);

        /* r_TE is negative, r_TM is positive and Lambda(l1,l2,m) is negative.
           => TM negative sign, TE positive sign */
        cint->signA_TM = -pow(-1, l2+m) * cint->signA_TM;
        cint->signA_TE = +pow(-1, l2+m) * cint->signA_TE;


        /* C */
        prefactor = log_m-tau-ln_tau; /* m*exp(-tau)/tau */
        cint->lnC_TE = prefactor + logadd_ms(lnC_TE, signs_C, N, &cint->signC_TE);
        cint->lnC_TM = prefactor + logadd_ms(lnC_TM, signs_C, N, &cint->signC_TM);

        cint->signC_TM = -pow(-1, l2+m+1) * cint->signC_TM;
        cint->signC_TE = +pow(-1, l2+m+1) * cint->signC_TE;


        /* D */
        /* prefactor is identical to C */
        cint->lnD_TE = prefactor + logadd_ms(lnD_TE, signs_D, N, &cint->signD_TE);
        cint->lnD_TM = prefactor + logadd_ms(lnD_TM, signs_D, N, &cint->signD_TM);

        cint->signD_TM = -pow(-1, l2+m) * cint->signD_TM;
        cint->signD_TE = +pow(-1, l2+m) * cint->signD_TE;
    }
    else
    {
        cint->lnA_TM = cint->lnA_TE = -INFINITY;
        cint->signA_TM = cint->signA_TE = +1;

        cint->lnC_TM = cint->lnC_TE = -INFINITY;
        cint->signC_TM = cint->signC_TE = +1;

        cint->lnD_TM = cint->lnD_TE = -INFINITY;
        cint->signD_TM = cint->signD_TE = +1;
    }

    xfree(ln_ABCD);
    xfree(signs_ABCD);
}

/*
* Multiply the polynomials p1 and p2. The result is stored in pdest, which
* must have at least size len_p1+len_p2-1
*/
void inline polymult(edouble p1[], size_t len_p1, edouble p2[], size_t len_p2, edouble pdest[])
{
    size_t power,i;

    for(power = 0; power < len_p1+len_p2-1; power++)
    {
        const int min = power-len_p2+1;

        pdest[power] = 0;
        for(i = MAX(0,min); i <= MIN(power,len_p1-1); i++)
            pdest[power] += p1[i]*p2[power-i];
    }
}

/* Integrate the function f(x)*exp(-x) from 0 to inf
* f(x) is the polynomial of length len stored in p
* l1,l2,m are needed to calculate the prefactor Lambda(l1,l2,m)
*
* This function returns the logarithm of the integral. The sign will be stored
* in sign.
*/
double log_polyintegrate(edouble p[], size_t len, int l1, int l2, int m, int *sign)
{
    size_t i;
    int sign_lnLambda;
    edouble value = 0;
    double lnLambda = casimir_lnLambda(l1, l2, m, &sign_lnLambda);
    double lnfac_max = lnfac(len-1);

    assert(!isnan(lnLambda));
    assert(!isinf(lnLambda));

    for(i = 0; i < len; i++)
        value += expq(lnfac(i)-lnfac_max)*p[i];

    assert(!isnan(value));
    assert(!isinf(value));
    *sign = (double)copysignq(1, value) * sign_lnLambda;
    return lnLambda+lnfac_max+logq(fabsq(value));
}

void polym(edouble p[], int m, edouble xi)
{
    size_t k;

    for(k = 0; k < m-1; k++)
        p[k] = 0;

    for(k = 0; k < m; k++)
        p[2*m-2-k] = binom(m-1,k)*pow(2*xi, k);
}

void polyplm(edouble pl1[], edouble pl2[], int l1, int l2, int m, edouble xi)
{
    int k;
    double log2xi = log(2*xi);

    for(k = 0; k <= l1-m; k++)
        pl1[k] = expq(lngamma(1+k+m+l1)-lngamma(1+l1-k-m)-lngamma(1+k)-lngamma(1+k+m)-k*log2xi);

    for(k = 0; k <= l2-m; k++)
        pl2[k] = expq(lngamma(1+k+m+l2)-lngamma(1+l2-k-m)-lngamma(1+k)-lngamma(1+k+m)-k*log2xi);
}

void polydplm(edouble pl1[], edouble pl2[], int l1, int l2, int m, edouble xi)
{
    int k;
    double log2xi = log(2*xi);

    if(m == 0)
    {
        for(k = 1; k <= l1; k++)
            pl1[k-1] = expq(lngamma(1+k+l1)-lngamma(1+k)-lngamma(1+l1-k)-lngamma(k)-(k-1)*log2xi)/2;

        for(k = 1; k <= l2; k++)
            pl2[k-1] = expq(lngamma(1+k+l2)-lngamma(1+k)-lngamma(1+l2-k)-lngamma(k)-(k-1)*log2xi)/2;
    }
    else
    {
        for(k = 0; k <= l1-m+1; k++)
            pl1[k] = 0;
        for(k = 0; k <= l2-m+1; k++)
            pl2[k] = 0;

        pl1[l1+1-m] += (l1-m+1)*expq(lngamma(2*l1+3)-lngamma(l1+2)-lngamma(l1-m+2)-(l1+1-m)*log2xi);
        for(k = 0; k <= l1-m; k++)
        {
            edouble common = expq(lngamma(1+k+l1+m)-lngamma(1+k+m)-lngamma(1+k)-lngamma(1+l1-k-m)-k*log2xi);
            pl1[k] += common*(pow_2(m)+m*(k-l1-1)-2.0*k*l1-2*k)/(m-l1+k-1);
            pl1[k+1] -= common*(l1+1)/xi;
        }

        pl2[l2+1-m] += (l2-m+1)*expq(lngamma(2*l2+3)-lngamma(l2+2)-lngamma(l2-m+2)-(l2+1-m)*log2xi);
        for(k = 0; k <= l2-m; k++)
        {
            edouble common = expq(lngamma(1+k+l2+m)-lngamma(1+k+m)-lngamma(1+k)-lngamma(1+l2-k-m)-k*log2xi);
            pl2[k] += common*(pow_2(m)+m*(k-l2-1)-2.0*k*l2-2*k)/(m-l2+k-1);
            pl2[k+1] -= common*(l2+1.)/xi;
        }
    }
}

/*
* Returns the integrals A,B,C,D for l1,l2,m,xi and p=TE,TM
*/
void casimir_integrate_perf(casimir_integrals_t *cint, int l1, int l2, int m, double nT)
{
    double lnA, lnB, lnC, lnD;
    int signA, signB, signC, signD;
    double xi = 2*nT;
    edouble pdpl1m[l1-m+2];
    edouble pdpl2m[l2-m+2];

    polydplm(pdpl1m,pdpl2m,l1,l2,m,xi);

    if(m == 0)
    {
        edouble pm[3];
        edouble interim[l1+2];
        edouble result[l1+l2+1];

        lnA = lnC = lnD = -INFINITY;
        signA = signC = signD = 0;

        polym(pm, 2, xi);

        polymult(pm, 3, pdpl1m, l1, interim);
        polymult(interim, l1+2, pdpl2m, l2, result);

        lnB = -xi-3*logq(xi)+log_polyintegrate(result, l1+l2+1, l1,l2,m, &signB);
        signB *= pow(-1, l2+1);
    }
    else
    {
        edouble logprefactor = -xi+logq(xi)-2*m*logq(xi)-m*logq(4);
        edouble pm[2*m-1];
        edouble ppl1m[l1-m+1];
        edouble ppl2m[l2-m+1];

        edouble pmppl1m[(sizeof(pm)+sizeof(ppl1m))/sizeof(edouble)-1];
        edouble pmpdpl1m[(sizeof(pm)+sizeof(pdpl1m))/sizeof(edouble)-1];

        edouble pmppl1mppl2m[(sizeof(pmppl1m)+sizeof(ppl2m))/sizeof(edouble)-1];
        edouble pmpdpl1mpdpl2m[(sizeof(pmpdpl1m)+sizeof(pdpl2m))/sizeof(edouble)-1];
        edouble pmppl1mpdpl2m[(sizeof(pmppl1m)+sizeof(pdpl2m))/sizeof(edouble)-1];
        edouble pmpdpl1mppl2m[(sizeof(pmpdpl1m)+sizeof(ppl2m))/sizeof(edouble)-1];

        polym(pm, m,xi);
        polyplm(ppl1m,ppl2m,l1,l2,m,xi);

        polymult(pm, sizeof(pm)/sizeof(edouble), ppl1m, sizeof(ppl1m)/sizeof(edouble), pmppl1m);
        polymult(pm, sizeof(pm)/sizeof(edouble), pdpl1m, sizeof(pdpl1m)/sizeof(edouble), pmpdpl1m);

        polymult(pmppl1m, sizeof(pmppl1m)/sizeof(edouble), ppl2m, sizeof(ppl2m)/sizeof(edouble), pmppl1mppl2m);
        polymult(pmppl1m, sizeof(pmppl1m)/sizeof(edouble), pdpl2m, sizeof(pdpl2m)/sizeof(edouble), pmppl1mpdpl2m);
        polymult(pmpdpl1m, sizeof(pmpdpl1m)/sizeof(edouble), ppl2m, sizeof(ppl2m)/sizeof(edouble), pmpdpl1mppl2m);
        polymult(pmpdpl1m, sizeof(pmpdpl1m)/sizeof(edouble), pdpl2m, sizeof(pdpl2m)/sizeof(edouble), pmpdpl1mpdpl2m);

        lnA = 2*logq(m)+logprefactor+log_polyintegrate(pmppl1mppl2m, sizeof(pmppl1mppl2m)/sizeof(edouble), l1,l2,m,&signA);
        signA *= pow(-1,l2);

        lnB = logprefactor+log_polyintegrate(pmpdpl1mpdpl2m, sizeof(pmpdpl1mpdpl2m)/sizeof(edouble), l1,l2,m,&signB);
        signB *= pow(-1,l2+1);
        
        lnC = logq(m)+logprefactor+log_polyintegrate(pmppl1mpdpl2m, sizeof(pmppl1mpdpl2m)/sizeof(edouble), l1,l2,m,&signC);
        signC *= pow(-1,l2+1);
        
        lnD = logq(m)+logprefactor+log_polyintegrate(pmpdpl1mppl2m, sizeof(pmpdpl1mppl2m)/sizeof(edouble), l1,l2,m,&signD);
        signD *= pow(-1,l2);
    }

    cint->lnA_TM   = lnA;
    cint->signA_TM = signA;

    cint->lnA_TE   = lnA;
    cint->signA_TE = -signA;

    cint->lnB_TM   = lnB;
    cint->signB_TM = signB;

    cint->lnB_TE   = lnB;
    cint->signB_TE = -signB;

    cint->lnC_TM   = lnC;
    cint->signC_TM = signC;

    cint->lnC_TE   = lnC;
    cint->signC_TE = -signC;

    cint->lnD_TM   = lnD;
    cint->signD_TM = signD;

    cint->lnD_TE   = lnD;
    cint->signD_TE = -signD;
}
