#define _ISOC99_SOURCE
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include "plm_fast.h"

#define Nlm(l,m) ( sqrt((2.0*l+1)/(l*(l+1)) * exp(gsl_sf_lngamma(1+l-m)-gsl_sf_lngamma(1+l+m))) )

static void _plm_array(int lmax, int m, double x, double *plm)
{
    int l;

    if(m == 0)
        plm[0] = 1;
    else
        plm[0] = pow(-1,m)*gsl_sf_doublefact(2*m-1)*pow(gsl_pow_2(x)-1, m*0.5); // l=m,m=m

    if(lmax == m)
        return;

    plm[1] = plm[0]*x*(2*m+1); // l=m+1, m=m

    for(l = m+2; l <= lmax; l++)
        plm[l-m] = ((2*l-1)*x*plm[l-m-1] - (l+m-1)*plm[l-m-2])/(l-m);
}

void YYY(int l1, int l2, int m, double x, double *yl1myl2m, double *dyl1mdyl2m, double *yl1mdyl2m, double *dyl1myl2m)
{
    int lmax = MAX(l1,l2)+1;
    double plm[lmax-m+1];
    double lambda_l1 = Nlm(l1,m);
    double lambda_l2;
    double denom = gsl_pow_2(x)-1;
    double Pl1m, Pl2m, dPl1m, dPl2m;

    if(l1 == l2)
        lambda_l2 = lambda_l1;
    else
        lambda_l2 = Nlm(l2,m);

    _plm_array(lmax, m, x, plm);

    // Lambda * Pl1m*Pl2m
    Pl1m = lambda_l1*plm[l1-m];
    Pl2m = lambda_l2*plm[l2-m];
    *yl1myl2m = pow(-1, 1+l2+m+(m%2))*Pl1m*Pl2m;

    // Lambda * dPl1m*dPl2m
    dPl1m = ((l1-m+1)*(lambda_l1*plm[l1+1-m]) - (l1+1)*x*(lambda_l1*plm[l1-m]))/denom;
    if(l1 == l2)
        dPl2m = dPl1m;
    else
        dPl2m = ((l2-m+1)*(lambda_l2*plm[l2+1-m]) - (l2+1)*x*(lambda_l2*plm[l2-m]))/denom;

    *dyl1mdyl2m = pow(-1, l2+m-(m%2))*dPl1m*dPl2m;
    
    *yl1mdyl2m = pow(-1, l2+m+1-(m%2))*Pl1m*dPl2m;
    *dyl1myl2m = pow(-1, l1+m+1-(m%2))*Pl2m*dPl1m;
}

double Yl1mYl2m(int l1, int l2, int m, double x)
{
    int lmax = MAX(l1,l2);
    double plm[lmax-m+1];
    double lambda_1 = Nlm(l1,m);
    double lambda_2;

    if(l1 == l2)
        lambda_2 = lambda_1;
    else
        lambda_2 = Nlm(l2,m);

    _plm_array(lmax, m, x, plm);

    return pow(-1, 1+l2+m+(m%2)) * lambda_1*plm[l1-m] * lambda_2*plm[l2-m];
}

double dYl1mdYl2m(int l1, int l2, int m, double x)
{
    int lmax = MAX(l1,l2)+1;
    double dPl1m, dPl2m;
    double denom = gsl_pow_2(x)-1;
    double plm[lmax-m+1];
    double lambda_l1 = Nlm(l1,m);

    _plm_array(lmax, m, x, plm);

    dPl1m = ((l1-m+1)*(lambda_l1*plm[l1+1-m]) - (l1+1)*x*(lambda_l1*plm[l1-m]))/denom;

    if(l1 == l2)
        dPl2m = dPl1m;
    else
    {
        double lambda_l2 = Nlm(l2,m);

        dPl2m = ((l2-m+1)*(lambda_l2*plm[l2+1-m]) - (l2+1)*x*(lambda_l2*plm[l2-m]))/denom;
    }

    return pow(-1, l2+m-(m%2))*dPl1m*dPl2m;
}

double Yl1mdYl2m(int l1, int l2, int m, double x)
{
    int lmax = MAX(l1,l2)+1;
    double Pl1m, dPl2m;
    double plm[lmax-m+1];
    double lambda_l1 = Nlm(l1,m);
    double lambda_l2;
    
    if(l1 == l2)
        lambda_l2 = lambda_l1;
    else
        lambda_l2 = Nlm(l2,m);

    _plm_array(lmax, m, x, plm);

    Pl1m = plm[l1-m]*lambda_l1;

    plm[l2+1-m] *= lambda_l2;
    plm[l2-m]   *= lambda_l2;

    dPl2m = ((l2-m+1)*plm[l2+1-m] - (l2+1)*x*plm[l2-m])/(gsl_pow_2(x)-1);

    return pow(-1, l2+m+1-(m%2))*Pl1m*dPl2m;
}
