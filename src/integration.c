#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "sfunc.h"
#include "libcasimir.h"
#include "integration.h"

#ifdef INTEGRATION_QUAD
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
        edouble value = 0;
        double lnLambda = casimir_lnLambda(l1, l2, m);
        double lnfac_max = lnfac(len-1);

        assert(!isnan(lnLambda));
        assert(!isinf(lnLambda));

        for(i = 0; i < len; i++)
            value += expq(lnfac(i)-lnfac_max)*p[i];

        assert(!isnan(value));
        assert(!isinf(value));
        *sign = (double)copysignq(1, value);
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
    void casimir_integrate(casimir_integrals_t *cint, int l1, int l2, int m, double xi)
    {
        edouble pdpl1m[l1-m+2];
        edouble pdpl2m[l2-m+2];

        polydplm(pdpl1m,pdpl2m,l1,l2,m,xi);

        if(m == 0)
        {
            edouble pm[3];
            edouble interim[l1+2];
            edouble result[l1+l2+1];

            cint->logA = cint->logC = cint->logD = -INFINITY;
            cint->signA = cint->signC = cint->signD = 0;

            polym(pm, 2, xi);

            polymult(pm, 3, pdpl1m, l1, interim);
            polymult(interim, l1+2, pdpl2m, l2, result);

            cint->logB = -xi-3*logq(xi)+log_polyintegrate(result, l1+l2+1, l1,l2,m,&cint->signB);
            cint->signB *= pow(-1, l2+1);
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

            cint->logA = 2*logq(m)+logprefactor+log_polyintegrate(pmppl1mppl2m, sizeof(pmppl1mppl2m)/sizeof(edouble), l1,l2,m,&cint->signA);
            cint->signA *= pow(-1,l2);

            cint->logB = logprefactor+log_polyintegrate(pmpdpl1mpdpl2m, sizeof(pmpdpl1mpdpl2m)/sizeof(edouble), l1,l2,m,&cint->signB);
            cint->signB *= pow(-1,l2+1);
            
            cint->logC = logq(m)+logprefactor+log_polyintegrate(pmppl1mpdpl2m, sizeof(pmppl1mpdpl2m)/sizeof(edouble), l1,l2,m,&cint->signC);
            cint->signC *= pow(-1,l2+1);
            
            cint->logD = logq(m)+logprefactor+log_polyintegrate(pmpdpl1mppl2m, sizeof(pmpdpl1mppl2m)/sizeof(edouble), l1,l2,m,&cint->signD);
            cint->signD *= pow(-1,l2);
        }
    }
#else
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
        size_t i, power;
        double list[len_p1+len_p2-1];
        char signs[len_p1+len_p2-1];

        for(power = 0; power < len_p1+len_p2-1; power++)
        {
            size_t len = 0;
            const int min = power-len_p2+1;
            for(i = MAX(0,min); i <= MIN(power,len_p1-1); i++)
            {
                const size_t j = power - i;

                if(p1[i].value != -INFINITY && p2[j].value != -INFINITY)
                {
                    list[len]  = p1[i].value + p2[j].value;
                    signs[len] = p1[i].sign  * p2[j].sign;
                    len++;
                }
            }

            if(len == 0)
            {
                pdest[power].value = -INFINITY;
                pdest[power].sign  = +1;
            }
            else
                pdest[power].value = logadd_ms(list, signs, len, &(pdest[power].sign));
        }
    }
#if 0
    void inline polymult(log_t p1[], size_t len_p1, log_t p2[], size_t len_p2, log_t pdest[])
    {
        size_t i, j;
        for(i = 0; i < len_p1+len_p2-1; i++)
        {
            pdest[i].value = -INFINITY;
            pdest[i].sign  = +1;
        }

        for(i = 0; i < len_p1; i++)
            for(j = 0; j < len_p2; j++)
                pdest[i+j].value = logadd_s(pdest[i+j].value, pdest[i+j].sign, p1[i].value+p2[j].value, p1[i].sign*p2[j].sign, &(pdest[i+j].sign));
    }
#endif

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
        double value;
        double lnLambda = casimir_lnLambda(l1, l2, m);

        value = p[0].value;
        *sign = p[0].sign;
        for(i = 1; i < len; i++)
            value = logadd_s(value, *sign, lnfac(i)+p[i].value, p[i].sign, sign);

        return lnLambda+value;
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
#endif
