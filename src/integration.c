#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "sfunc.h"
#include "libcasimir.h"
#include "integration.h"

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
void casimir_integrate(casimir_integrals_t *cint, int l1, int l2, int m, double nT)
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
