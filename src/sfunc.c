#include <math.h>
#include <quadmath.h>

#include "sfunc.h"


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
