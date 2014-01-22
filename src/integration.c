#include <complex.h>
#include <math.h>
#include <quadmath.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "sfunc.h"
#include "libcasimir.h"
#include "integration.h"

void polyprint(log_t p[], size_t len);

void polyprint(log_t p[], size_t len)
{
    int k;

    for(k = 0; k < len; k++)
        printf("%+gx^%d ", p[k].sign*exp(p[k].value), k);
    printf("\n");
}

/*
 * Multiply the polynomials p1 and p2. The result is stored in pdest, which
 * must have at least size len_p1+len_p2-1
 */
void inline polymult(log_t p1[], size_t len_p1, log_t p2[], size_t len_p2, log_t pdest[])
{
    size_t i, j, N, len;
    double *list_values = __builtin_assume_aligned(malloc(sizeof(double)*(len_p1+len_p2)), 16);
    int *list_signs     = malloc(sizeof(int)*(len_p1+len_p2));
    double max;

    for(N = 0; N < len_p1+len_p2-1; N++)
    {
        len = 0;
        max = -INFINITY;
        for(i = 0; i <= MIN(N,len_p1-1); i++)
        {
            j = N-i;
            if(j >= len_p2)
                continue;

            list_values[len] = p1[i].value+p2[j].value;
            list_signs[len]  = p1[i].sign*p2[j].sign;
            max = MAX(list_values[len], max);

            len++;
        }

        if(isinf(max))
        {
            pdest[N].value = -INFINITY;
            pdest[N].sign = 1;
        }
        else
        {
            double sum = 0;
            for(i = 0; i < len; i++)
            {
                sum += list_signs[i]*exp(list_values[i]-max);
                //printf("sum %d, %d %g => %g\n", (int)i, list[i].sign, list[i].value, sum);
            }

            pdest[N].value = max+log(fabs(sum));
            pdest[N].sign = copysign(1, sum);
        }
    }

    free(list_values);
    free(list_signs);
}

/* Integrate the function f(x)*exp(-x) from 0 to inf 
 * f(x) is the polynomial of length len stored in p
 * l1,l2,m are needed to calculate the prefactor Lambda(l1,l2,m)
 * 
 * This function returns the logarithm of the integral. The sign will be stored
 * in sign.
 */
double log_polyintegrate(log_t p[], size_t len, int l1, int l2, int m, int *sign)
{
    size_t i;
    double value = 0;
    double lnLambda = casimir_lnLambda(l1, l2, m);
    double lnfac_max = lnfac(len-1);

    assert(!isnan(lnLambda));
    assert(!isinf(lnLambda));

    for(i = 0; i < len; i++)
        value += exp(lnfac(i)+p[i].value-lnfac_max)* p[i].sign;

    *sign = copysign(1, value);

    assert(!isnan(value));
    //assert(!isinf(value));
    return lnLambda+lnfac_max+log(fabs(value));
}

void polym(log_t p[], int m, double xi)
{
    size_t k;
    double log2xi = log(2*xi);

    for(k = 0; k < m-1; k++)
    {
        p[k].value = -INFINITY;
        p[k].sign = +1;
    }

    for(k = 0; k < m; k++)
    {
        p[2*m-2-k].value = lbinom(m-1,k)+k*log2xi;
        p[2*m-2-k].sign = +1;
    }
}

void polyplm(log_t pl1[], log_t pl2[], int l1, int l2, int m, double xi)
{
    int k;
    double log2xi = log(2*xi);

    for(k = 0; k <= l1-m; k++)
    {
        pl1[k].value = lngamma(1+k+m+l1)-lngamma(1+l1-k-m)-lngamma(1+k)-lngamma(1+k+m)-k*log2xi;
        pl1[k].sign = +1;
    }

    for(k = 0; k <= l2-m; k++)
    {
        pl2[k].value = lngamma(1+k+m+l2)-lngamma(1+l2-k-m)-lngamma(1+k)-lngamma(1+k+m)-k*log2xi;
        pl2[k].sign = +1;
    }
}

void polydplm(log_t pl1[], log_t pl2[], int l1, int l2, int m, double xi)
{
    int k;
    double log2xi = log(2*xi);

    if(m == 0)
    {
        for(k = 1; k <= l1; k++)
        {
            pl1[k-1].value = lngamma(1+k+l1)-lngamma(1+k)-lngamma(1+l1-k)-lngamma(k)-(k-1)*log2xi-M_LN2;
            pl1[k-1].sign = +1;
        }

        for(k = 1; k <= l2; k++)
        {
            pl2[k-1].value = lngamma(1+k+l2)-lngamma(1+k)-lngamma(1+l2-k)-lngamma(k)-(k-1)*log2xi-M_LN2;
            pl2[k-1].sign = +1;
        }
    }
    else
    {
        for(k = 0; k <= l1-m+1; k++)
        {
            pl1[k].value= -INFINITY;
            pl1[k].sign= +1;
        }
        for(k = 0; k <= l2-m+1; k++)
        {
            pl2[k].value= -INFINITY;
            pl2[k].sign= +1;
        }

        //pl1[l1+1-m] += (l1-m+1)*exp(lngamma(2*l1+3)-lngamma(l1+2)-lngamma(l1-m+2)-(l1+1-m)*log2xi);
        pl1[l1+1-m].value = logadd_s(pl1[l1+1-m].value, pl1[l1+1-m].sign, log(l1-m+1)+lngamma(2*l1+3)-lngamma(l1+2)-lngamma(l1-m+2)-(l1+1-m)*log2xi, 1, &(pl1[l1+1-m].sign));
        for(k = 0; k <= l1-m; k++)
        {
            double common = lngamma(1+k+l1+m)-lngamma(1+k+m)-lngamma(1+k)-lngamma(1+l1-k-m)-k*log2xi;
            double fac = (pow_2(m)+m*(k-l1-1)-2.*k*l1-2*k)/(m-l1+k-1);
            //pl1[k]   += common*(pow_2(m)+m*(k-l1-1)-2*k*l1-2*k)/(m-l1+k-1);
            pl1[k].value = logadd_s(pl1[k].value, pl1[k].sign, common+log(fabs(fac)), copysign(1,fac), &(pl1[k].sign));
            //pl1[k+1] -= common*(l1+1)/xi;
            pl1[k+1].value = logadd_s(pl1[k+1].value, pl1[k+1].sign, common+log(l1+1)-log(xi), -1, &(pl1[k+1].sign));
        }

        //pl2[l2+1-m] += (l2-m+1)*exp(lngamma(2*l2+3)-lngamma(l2+2)-lngamma(l2-m+2)-(l2+1-m)*log2xi);
        pl2[l2+1-m].value = logadd_s(pl2[l2+1-m].value, pl2[l2+1-m].sign, log(l2-m+1)+lngamma(2*l2+3)-lngamma(l2+2)-lngamma(l2-m+2)-(l2+1-m)*log2xi, 1, &(pl2[l2+1-m].sign));
        for(k = 0; k <= l2-m; k++)
        {
            double common = lngamma(1+k+l2+m)-lngamma(1+k+m)-lngamma(1+k)-lngamma(1+l2-k-m)-k*log2xi;
            double fac = (pow_2(m)+m*(k-l2-1)-2.*k*l2-2*k)/(m-l2+k-1);

            //pl2[k]   += common*(pow_2(m)+m*(k-l2-1)-2.*k*l2-2*k)/(m-l2+k-1);
            pl2[k].value = logadd_s(pl2[k].value, pl2[k].sign, common+log(fabs(fac)), copysign(1,fac), &(pl2[k].sign));

            //pl2[k+1] -= common*(l2+1.)/xi;
            pl2[k+1].value = logadd_s(pl2[k+1].value, pl2[k+1].sign, common+log(l2+1)-log(xi), -1, &(pl2[k+1].sign));
        }
    }
}

/*
 * Returns the integrals A,B,C,D for l1,l2,m,xi and p=TE,TM
 */
void casimir_integrate(casimir_integrals_t *cint, int l1, int l2, int m, double xi)
{
    log_t pdpl1m[l1-m+2];
    log_t pdpl2m[l2-m+2];

    polydplm(pdpl1m,pdpl2m,l1,l2,m,xi);

    if(m == 0)
    {
        log_t pm[3];
        log_t interim[l1+2];
        log_t result[l1+l2+1];

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
        log_t pm[2*m-1];
        log_t ppl1m[l1-m+1];
        log_t ppl2m[l2-m+1];

        log_t pmppl1m[(sizeof(pm)+sizeof(ppl1m))/sizeof(log_t)-1];
        log_t pmpdpl1m[(sizeof(pm)+sizeof(pdpl1m))/sizeof(log_t)-1];

        log_t pmppl1mppl2m[(sizeof(pmppl1m)+sizeof(ppl2m))/sizeof(log_t)-1];
        log_t pmpdpl1mpdpl2m[(sizeof(pmpdpl1m)+sizeof(pdpl2m))/sizeof(log_t)-1];
        log_t pmppl1mpdpl2m[(sizeof(pmppl1m)+sizeof(pdpl2m))/sizeof(log_t)-1];
        log_t pmpdpl1mppl2m[(sizeof(pmpdpl1m)+sizeof(ppl2m))/sizeof(log_t)-1];

        polym(pm, m,xi);
        polyplm(ppl1m,ppl2m,l1,l2,m,xi);

        polymult(pm, sizeof(pm)/sizeof(log_t), ppl1m,  sizeof(ppl1m)/sizeof(log_t), pmppl1m);
        polymult(pm, sizeof(pm)/sizeof(log_t), pdpl1m, sizeof(pdpl1m)/sizeof(log_t), pmpdpl1m);

        polymult(pmppl1m,  sizeof(pmppl1m)/sizeof(log_t),  ppl2m,  sizeof(ppl2m)/sizeof(log_t),  pmppl1mppl2m);
        polymult(pmppl1m,  sizeof(pmppl1m)/sizeof(log_t),  pdpl2m, sizeof(pdpl2m)/sizeof(log_t), pmppl1mpdpl2m);
        polymult(pmpdpl1m, sizeof(pmpdpl1m)/sizeof(log_t), ppl2m,  sizeof(ppl2m)/sizeof(log_t),  pmpdpl1mppl2m);
        polymult(pmpdpl1m, sizeof(pmpdpl1m)/sizeof(log_t), pdpl2m, sizeof(pdpl2m)/sizeof(log_t), pmpdpl1mpdpl2m);

        cint->logA = 2*log(m)+logprefactor+log_polyintegrate(pmppl1mppl2m, sizeof(pmppl1mppl2m)/sizeof(log_t), l1,l2,m,&cint->signA);
        cint->signA *= pow(-1,l2);

        cint->logB = logprefactor+log_polyintegrate(pmpdpl1mpdpl2m, sizeof(pmpdpl1mpdpl2m)/sizeof(log_t), l1,l2,m,&cint->signB);
        cint->signB *= pow(-1,l2+1);
        
        cint->logC = log(m)+logprefactor+log_polyintegrate(pmppl1mpdpl2m, sizeof(pmppl1mpdpl2m)/sizeof(log_t), l1,l2,m,&cint->signC);
        cint->signC *= pow(-1,l2+1);
        
        cint->logD = log(m)+logprefactor+log_polyintegrate(pmpdpl1mppl2m, sizeof(pmpdpl1mppl2m)/sizeof(log_t), l1,l2,m,&cint->signD);
        cint->signD *= pow(-1,l2);
    }
}
