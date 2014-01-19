#include <math.h>
#include <quadmath.h>

#include "sfunc.h"


double bessel_lnInuKnu(const int nu, const double x, double *lnInu_p, double *lnKnu_p)
{
    int k;
    __float128 lnInu, lnKnu, lnKnup;

    // calculate Knu and Knu+
    {
        double prefactor = -x+0.5*(M_LOGPI-M_LN2-log(x));
        __float128 Km = 1, Kp = 1+1./x;

        if(nu == 0)
            return prefactor+logq(Km);
        else if(nu == 1)
            return prefactor+logq(Kp);
        else
            for(k = 2; k <= nu+1; k++)
            {
                __float128 Kn = (2*k-1)*Kp/x + Km;
                Km = Kp;
                Kp = Kn;
            }

        lnKnu  = prefactor+logq(Km);
        lnKnup = prefactor+logq(Kp);
    }

    {
        #define an(n,nu,x) (2*(nu+0.5+n)/x)

        double nom   = an(2,nu,x)+1/an(1,nu,x);
        double denom = an(2,nu,x);
        double ratio = (an(1,nu,x)*nom)/denom;
        double ratio_last = 0;

        k = 3;
        while(1)
        {
            nom   = an(k,nu,x)+1/nom;
            denom = an(k,nu,x)+1/denom;
            ratio *= nom/denom;

            if(ratio_last != 0 && fabs(1-ratio/ratio_last) < 1e-10)
                break;

            ratio_last = ratio;
            k++;
        }

        return -log(x*(exp(lnKnup)+1/ratio*exp(lnKnu)));
    }

    if(lnInu_p != NULL)
        *lnInu_p = lnInu;
    if(lnKnu_p != NULL)
        *lnKnu_p = lnKnu;
}

#if 1
double bessel_lnInu(const int n, const double x)
{
    int k = 3;
    const double nu = n+0.5;

    #undef an
    #define an(n,nu,x) (2*(nu+n)/x)

    double nom   = an(2,nu,x)+1/an(1,nu,x);
    double denom = an(2,nu,x);
    double ratio = (an(1,nu,x)*nom)/denom;
    double ratio_last = 0;
    const double lnKnu  = bessel_lnKnu(n,x);
    const double lnKnup = bessel_lnKnu(n+1,x);

    while(1)
    {
        nom   = an(k,nu,x)+1/nom;
        denom = an(k,nu,x)+1/denom;
        ratio *= nom/denom;

        if(ratio_last != 0 && fabs(1-ratio/ratio_last) < 1e-10)
            break;

        ratio_last = ratio;
        k++;
    }

    return -log(x*(exp(lnKnup)+1/ratio*exp(lnKnu)));
}
#else
double bessel_lnInu(const int n, const double x)
{
    double nu = n+0.5;

    if(x < sqrt(nu+1)/1e4)
        return nu*(log(x)-M_LN2)-lgamma(nu+1);
    else if(x > nu*nu*1e4)
        return x-0.5*log(2*M_PI*x)+log1p((4*nu*nu-1)/(8*x));
    else
    {
        int k = 0;
        const double prefactor = nu*(log(x)-M_LN2);
        const __float128 xsqr4 = 2*(logq(x)-M_LN2q);
        __float128 value = 0;

        while(1)
        {
            __float128 denom = lgammaq(1+k)+lgammaq(nu+k+1);
            __float128 term = expq(k*xsqr4-denom);
            value += term;

            if(fabsq(term/value) < 1e-25)
                break;

            k++;
        }

        return prefactor+logq(fabsq(value));
    }
}
#endif

double bessel_lnKnu(const int nu, const double x)
{
    int l;
    double prefactor = -x+0.5*(M_LOGPI-M_LN2-log(x));
    __float128 Km = 1;
    __float128 Kp = 1+1./x;

    if(nu == 0)
        return prefactor+logq(Km);
    else if(nu == 1)
        return prefactor+logq(Kp);

    for(l = 2; l <= nu; l++)
    {
        __float128 Kn = (2*l-1)*Kp/x + Km;
        Km = Kp;
        Kp = Kn;
    }

    return prefactor+logq(Kp);
}
