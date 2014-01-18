#include <math.h>
#include <quadmath.h>

#include "sfunc.h"

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

            if(fabsq(term/value) < 1e-20)
                break;

            k++;
        }

        return prefactor+logq(fabsq(value));
    }
}

double bessel_lnKnu(const int nu, const double x)
{
    int l;
    double prefactor = -x+0.5*(M_LOGPIE-M_LN2-log(x));
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
