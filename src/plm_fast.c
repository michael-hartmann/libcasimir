#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include "plm_fast.h"

/* This module implements associated legendre functions and its derivatives
 * for m >= 0 and x >= 1.
 * 
 * Associated Legendre polynomials are defined as follows:
 *     Plm(x) = (-1)^m (1-x²)^(m/2) * d^m/dx^m Pl(x)
 * where Pl(x) denotes a Legendre polynomial.
 *
 * As Pl(x) are ordinary polynomials, the only problem is the term (1-x²) when
 * extending the domain to values of x > 1. The continuation is ambiguous.
 * We will implement 
 *     Plm(x) = (-1)^m (x²-1)^(m/2) * d^m/dx^m Pl(x)
 * here.
 *
 * Spherical harmonics at phi=0 are just associated Legendre functions
 * multiplied by a constant:
 *     Ylm(x,0) = Nlm * Plm(cos(theta))
 * where
 *     Nlm = sqrt( (2*l+1) * (l-m)!/(l+m)! ).
 * The functions assume x=cos(x) and calculate
 *     Ylm(x,0) = Nlm * Plm(x).
 *
 * (*) Note:
 * Products of associated legendre polynomials with common m are unambiguous, because
 *     (i)² = (-i)² = -1.
 */

/* calculate Plm for l=m...l=lmax */
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

/* calculate Plm(x) */
double plm_Plm(int l, int m, double x)
{
    double plm[l-m+1];
    _plm_array(l, m, x, plm);
    return plm[l-m];
}

/* calculate dPlm(x) */
double plm_dPlm(int l, int m, double x)
{
    int lmax = l+1;
    double denom = gsl_pow_2(x)-1;
    double plm[lmax-m+1];

    _plm_array(lmax, m, x, plm);

    return ((l-m+1)*(Nlm(l,m)*plm[l+1-m]) - (l+1)*x*(Nlm(l,m)*plm[l-m]))/denom;
}

/* calculate Pl1m(x), Pl2m(x), dPl1m(x), dPl2m(x) */
void plm_Yl12md(int l1, int l2, int m, double x, double *Yl1m, double *Yl2m, double *dYl1m, double *dYl2m)
{
    int lmax = MAX(l1,l2)+1;
    double plm[lmax-m+1];
    double lambda_l1 = Nlm(l1,m);
    double lambda_l2;
    double denom = gsl_pow_2(x)-1;

    if(l1 == l2)
        lambda_l2 = lambda_l1;
    else
        lambda_l2 = Nlm(l2,m);

    _plm_array(lmax, m, x, plm);

    *Yl1m = lambda_l1*plm[l1-m];
    *Yl2m = lambda_l2*plm[l2-m];

    *dYl1m = ((l1-m+1)*(lambda_l1*plm[l1+1-m]) - (l1+1)*x*(lambda_l1*plm[l1-m]))/denom;
    if(l1 == l2)
        *dYl2m = *dYl1m;
    else
        *dYl2m = ((l2-m+1)*(lambda_l2*plm[l2+1-m]) - (l2+1)*x*(lambda_l2*plm[l2-m]))/denom;
}

/* calcualte Yl1m(x)*Yl2m(x), see (*) */
double plm_Yl1mYl2m(int l1, int l2, int m, double x)
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

/* calculate dYl1m(x)*Yl2m(x), see (*) */
double plm_dYl1mdYl2m(int l1, int l2, int m, double x)
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

/* calculate Yl1m(x)*Yl2m(x), see (*) */
double plm_Yl1mdYl2m(int l1, int l2, int m, double x)
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
