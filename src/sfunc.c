#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "quad.h"
#include "sfunc.h"

double inline logadd_s(const double a, const int sign_a, const double b, const int sign_b, int *sign)
{
    if(a == -INFINITY)
    {
        *sign = sign_b;
        return b;
    }
    else if(b == -INFINITY)
    {
        *sign = sign_a;
        return a;
    }

    if(a > b)
    {
        *sign = sign_a;
        return a + log1p(sign_a*sign_b*exp(b-a));
    }
    else
    {
        *sign = sign_b;
        return b + log1p(sign_a*sign_b*exp(a-b));
    }
}

double inline logadd(const double a, const double b)
{
    if(a == -INFINITY)
        return b;
    else if(b == -INFINITY)
        return a;
    else if(a < b)
        return b + log1p(exp(a-b));
    else
        return a + log1p(exp(b-a));
}

double inline logadd_m(const double list[], size_t len)
{
    size_t i;
    double sum;
    double max = list[0];

    for(i = 1; i < len; i++)
        if(list[i] > max)
            max = list[i];

    sum = exp(list[0]-max);
    for(i = 1; i < len; i++)
        sum += exp(list[i]-max);

    return max + log(sum);
}

double inline logadd_ms(const double list[], const char signs[], const size_t len, int *sign)
{
    size_t i;
    double sum;
    double max = list[0];

    for(i = 1; i < len; i++)
        if(list[i] > max)
            max = list[i];

    sum = signs[0]*exp(list[0]-max);
    for(i = 1; i < len; i++)
        sum += signs[i]*exp(list[i]-max);

    *sign = copysign(1, sum);
    return max + log(fabs(sum));
}

double inline lbinom(int n, int k)
{
    return lngamma(1+n)-lngamma(1+k)-lngamma(1+n-k);
}

double inline binom(int n, int k)
{
    return exp(lngamma(1+n)-lngamma(1+k)-lngamma(1+n-k));
}

void bessel_lnInuKnu(int nu, const double x, double *lnInu_p, double *lnKnu_p)
{
    int l;
    __float128 lnKnu = 1, lnKnup = 1+1./x;

    // calculate Knu, Knup
    {
        double prefactor = -x+0.5*(M_LOGPI-M_LN2-log(x));

        if(nu == 0)
        {
            lnKnu  = prefactor+logq(lnKnu);
            lnKnup = prefactor+logq(lnKnup);
        }
        else
        {
            for(l = 2; l <= nu+1; l++)
            {
                __float128 Kn = (2*l-1)*lnKnup/x + lnKnu;
                lnKnu  = lnKnup;
                lnKnup = Kn;
            }

            lnKnup = prefactor+logq(lnKnup);
            lnKnu  = prefactor+logq(lnKnu);
        }

        if(isnanq(lnKnup) || isinfq(lnKnup))
        {
            /* so, we couldn't calculate lnKnup and lnKnu. Maybe we can at
             * least use the asymptotic behaviour for small values.
             */
            if(x < sqrt(nu)*1e3)
            {
                /* small arguments */
                lnKnu  = lngamma(nu+0.5)-M_LN2+(nu+0.5)*(M_LN2-log(x));
                lnKnup = lngamma(nu+1.5)-M_LN2+(nu+1.5)*(M_LN2-log(x));
            }
            else
                lnKnu = lnKnup = 0;
        }

        if(lnKnu_p != NULL)
            *lnKnu_p = lnKnu;
    }

    if(lnInu_p != NULL)
    {
        #define an(n,nu,x) (2*(nu+0.5+n)/x)

        __float128 nom   = an(2,nu,x)+1/an(1,nu,x);
        __float128 denom = an(2,nu,x);
        __float128 ratio = (an(1,nu,x)*nom)/denom;
        __float128 ratio_last = 0;

        l = 3;
        while(1)
        {
            nom   = an(l,nu,x)+1/nom;
            denom = an(l,nu,x)+1/denom;
            ratio *= nom/denom;

            if(ratio_last != 0 && fabs(1-ratio/ratio_last) < 1e-10)
                break;

            ratio_last = ratio;
            l++;
        }

        *lnInu_p = -log(x)-lnKnu-logq(expq(lnKnup-lnKnu)+1/ratio);
    }
}

double bessel_lnKnu(const int nu, const double x)
{
    double Knu;
    bessel_lnInuKnu(nu, x, NULL, &Knu);
    return Knu;
}

double bessel_lnInu(const int nu, const double x)
{
    double Inu;
    bessel_lnInuKnu(nu, x, &Inu, NULL);
    return Inu;
}
