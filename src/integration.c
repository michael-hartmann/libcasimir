#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "sfunc.h"
#include "libcasimir.h"
#include "integration.h"

void polyprint(double p[], size_t len);

void polyprint(double p[], size_t len)
{
    int k;

    for(k = 0; k < len; k++)
        if(p[k] != 0)
            printf("%+gx^%d ", p[k], k);
    printf("\n");
}

/*
 * Multiply the polynomials p1 and p2. The result is stored in pdest, which
 * must have at least size len_p1+len_p2-1
 */
void inline polymult(double p1[], size_t len_p1, double p2[], size_t len_p2, double pdest[])
{
    size_t i, j;
    for(i = 0; i < len_p1+len_p2-1; i++)
        pdest[i] = 0;

    for(i = 0; i < len_p1; i++)
        for(j = 0; j < len_p2; j++)
            pdest[i+j] += p1[i]*p2[j];
}

double polyintegrate(double p[], size_t len, int l1, int l2, int m, double scale)
{
    size_t i;
    double value = 0;
    double lnLambda = casimir_lnLambda(l1, l2, m);
    double logscale = log(fabs(scale));
    double sign     = copysign(1, scale);

    for(i = 0; i < len; i++)
        value += sign*exp(logscale+lnLambda+lngamma(1+i))*p[i];

    return value;
}

double log_polyintegrate(double p[], size_t len, int l1, int l2, int m, int *sign)
{
    size_t i;
    double value = 0;
    double lnLambda = casimir_lnLambda(l1, l2, m);
    double lnfac_max = lnfac(len-1);

    for(i = 0; i < len; i++)
        value += exp(lnfac(i)-lnfac_max)*p[i];

    *sign = copysign(1, value);
    return lnLambda+lnfac_max+log(fabs(value));
}

void polym(double p[], int m, double xi)
{
    size_t k;

    for(k = 0; k < m-1; k++)
        p[k] = 0;

    for(k = 0; k < m; k++)
        p[2*m-2-k] = binom(m-1,k)*pow(2*xi, k);
}

void polyplm(double pl1[], double pl2[], int l1, int l2, int m, double xi)
{
    int k;
    double log2xi = log(2*xi);

    for(k = 0; k <= l1-m; k++)
        pl1[k] = exp(lngamma(1+k+m+l1)-lngamma(1+l1-k-m)-lngamma(1+k)-lngamma(1+k+m)-k*log2xi);

    for(k = 0; k <= l2-m; k++)
        pl2[k] = exp(lngamma(1+k+m+l2)-lngamma(1+l2-k-m)-lngamma(1+k)-lngamma(1+k+m)-k*log2xi);
}

void polydplm(double pl1[], double pl2[], int l1, int l2, int m, double xi)
{
    int k;
    double log2xi = log(2*xi);

    if(m == 0)
    {
        for(k = 1; k <= l1; k++)
            pl1[k-1] = exp(lngamma(1+k+l1)-lngamma(1+k)-lngamma(1+l1-k)-lngamma(k)-(k-1)*log2xi)/2;

        for(k = 1; k <= l2; k++)
            pl2[k-1] = exp(lngamma(1+k+l2)-lngamma(1+k)-lngamma(1+l2-k)-lngamma(k)-(k-1)*log2xi)/2;
    }
    else
    {
        for(k = 0; k <= l1-m+1; k++)
            pl1[k] = 0;
        for(k = 0; k <= l2-m+1; k++)
            pl2[k] = 0;

        pl1[l1+1-m] += (l1-m+1)*exp(lngamma(2*l1+3)-lngamma(l1+2)-lngamma(l1-m+2)-(l1+1-m)*log2xi);
        for(k = 0; k <= l1-m; k++)
        {
            double common = exp(lngamma(1+k+l1+m)-lngamma(1+k+m)-lngamma(1+k)-lngamma(1+l1-k-m)-k*log2xi);
            pl1[k]   += common*(pow_2(m)+m*(k-l1-1)-2*k*l1-2*k)/(m-l1+k-1);
            pl1[k+1] -= common*(l1+1)/xi;
        }

        pl2[l2+1-m] += (l2-m+1)*exp(lngamma(2*l2+3)-lngamma(l2+2)-lngamma(l2-m+2)-(l2+1-m)*log2xi);
        for(k = 0; k <= l2-m; k++)
        {
            double common = exp(lngamma(1+k+l2+m)-lngamma(1+k+m)-lngamma(1+k)-lngamma(1+l2-k-m)-k*log2xi);
            pl2[k]   += common*(pow_2(m)+m*(k-l2-1)-2.0*k*l2-2*k)/(m-l2+k-1);
            pl2[k+1] -= common*(l2+1.)/xi;
        }
    }
}

/*
 * Returns the integrals A,B,C,D for l1,l2,m,xi and p=TE,TM
 */
int casimir_integrate(casimir_integrals_t *cint, int l1, int l2, int m, double xi)
{
    double pdpl1m[l1-m+2];
    double pdpl2m[l2-m+2];

    polydplm(pdpl1m,pdpl2m,l1,l2,m,xi);

    if(m == 0)
    {
        double pm[3];
        double interim[l1+2];
        double result[l1+l2+1];

        cint->logA  = cint->logC  = cint->logD  = -INFINITY;
        cint->signA = cint->signC = cint->signD = 0;

        polym(pm, 2, xi);

        polymult(pm, 3, pdpl1m, l1, interim);
        polymult(interim, l1+2, pdpl2m, l2, result);

        cint->logB = -xi-3*log(xi)+log_polyintegrate(result, l1+l2+1, l1,l2,m,&cint->signB);
        cint->signB *= pow(-1, l2+1);
    }
    else
    {
        double logprefactor = -xi+log(xi)-2*m*log(xi)-m*log(4);
        double pm[2*m-1];
        double ppl1m[l1-m+1];
        double ppl2m[l2-m+1];

        double pmppl1m[(sizeof(pm)+sizeof(ppl1m))/sizeof(double)-1];
        double pmpdpl1m[(sizeof(pm)+sizeof(pdpl1m))/sizeof(double)-1];

        double pmppl1mppl2m[(sizeof(pmppl1m)+sizeof(ppl2m))/sizeof(double)-1];
        double pmpdpl1mpdpl2m[(sizeof(pmpdpl1m)+sizeof(pdpl2m))/sizeof(double)-1];
        double pmppl1mpdpl2m[(sizeof(pmppl1m)+sizeof(pdpl2m))/sizeof(double)-1];
        double pmpdpl1mppl2m[(sizeof(pmpdpl1m)+sizeof(ppl2m))/sizeof(double)-1];

        polym(pm, m,xi);
        polyplm(ppl1m,ppl2m,l1,l2,m,xi);

        polymult(pm, sizeof(pm)/sizeof(double), ppl1m,  sizeof(ppl1m)/sizeof(double), pmppl1m);
        polymult(pm, sizeof(pm)/sizeof(double), pdpl1m, sizeof(pdpl1m)/sizeof(double), pmpdpl1m);

        polymult(pmppl1m,  sizeof(pmppl1m)/sizeof(double),  ppl2m,  sizeof(ppl2m)/sizeof(double),  pmppl1mppl2m);
        polymult(pmppl1m,  sizeof(pmppl1m)/sizeof(double),  pdpl2m, sizeof(pdpl2m)/sizeof(double), pmppl1mpdpl2m);
        polymult(pmpdpl1m, sizeof(pmpdpl1m)/sizeof(double), ppl2m,  sizeof(ppl2m)/sizeof(double),  pmpdpl1mppl2m);
        polymult(pmpdpl1m, sizeof(pmpdpl1m)/sizeof(double), pdpl2m, sizeof(pdpl2m)/sizeof(double), pmpdpl1mpdpl2m);

        cint->logA = 2*log(m)+logprefactor+log_polyintegrate(pmppl1mppl2m,   sizeof(pmppl1mppl2m)/sizeof(double),   l1,l2,m,&cint->signA);
        cint->signA *= pow(-1,l2);

        cint->logB = logprefactor+log_polyintegrate(pmpdpl1mpdpl2m, sizeof(pmpdpl1mpdpl2m)/sizeof(double), l1,l2,m,&cint->signB);
        cint->signB *= pow(-1,l2+1);
        
        cint->logC = log(m)+logprefactor+log_polyintegrate(pmppl1mpdpl2m,  sizeof(pmppl1mpdpl2m)/sizeof(double),  l1,l2,m,&cint->signC);
        cint->signC *= pow(-1,l2+1);
        
        cint->logD = log(m)+logprefactor+log_polyintegrate(pmpdpl1mppl2m,  sizeof(pmpdpl1mppl2m)/sizeof(double),  l1,l2,m,&cint->signD);
        cint->signD *= pow(-1,l2);
    }

    return 0;
}
