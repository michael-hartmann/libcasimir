#define _ISOC99_SOURCE
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>

#include "plm.h"

double Nlm(int l, int m)
{
    return sqrt((2*l+1) * exp(gsl_sf_lngamma(1+l-m)-gsl_sf_lngamma(1+l+m)));
}

/* ============================================================================
 * Purpose: Compute the associated Legendre functions Plm(x) and their
 * derivatives Plm'(z) for a real argument argument
 *
 * Input:
 *     x:     argument of associated Legendre function
 *     l:     degree of Plm(x), l = 0,1,2,...,n
 *     m:     order of Plm(x), m = 0,1,2,...,l
 *     pmPlm:  
 *     pmdPlm: 
 *
 * Return:
 *     value of Plm(x)
 *
 * Restrictions:
 *     m >= 0
 *     l >= m
 *     x >= 0
 * https://github.com/scipy/scipy/blob/master/scipy/special/specfun/specfun.f
 * ============================================================================
 */
double Plm(int l, int m, double x, gsl_matrix **pmPlm, gsl_matrix **pmdPlm)
{
    int il,im;
    double xr,xs,result;
    gsl_matrix *mPlm;

    mPlm = gsl_matrix_alloc(l+1, l+1);
    gsl_matrix_set_zero(mPlm);
    if(pmdPlm != NULL)
    {
        *pmdPlm = gsl_matrix_alloc(l+1, l+1);
        gsl_matrix_set_zero(*pmdPlm);
    }

    /* Plm(0,0,x) = 1 */
    gsl_matrix_set(mPlm,0,0,1);

    /* if l == 0, we're done */
    if(l == 0)
        return 1;

    /* if x == 1, things are different */
    if(fabs(x) == 1)
    {
        /* Plm(l,0,x) = 1 */
        for(il = 1; il <= l; il++)
            gsl_matrix_set(mPlm, il, 0, 1);

        if(pmdPlm != NULL)
        {
            /* Plm'(l,0,x) = m*(m+1)/2 */
            for(il = 1; il <= l; il++)
                gsl_matrix_set(*pmdPlm, il, 0, m*(m+1)/2);

            if(m >= 1)
            {
                for(il = 1; il <= l; il++)
                {
                    /* Plm'(l,1,x) = m*(m+1)/2 */
                    gsl_matrix_set(*pmdPlm, il, 1, INFINITY);
                    if(m >= 2)
                        /* Plm'(l,2,x) = -(l+2)(l+1)*l*(l-1)/4 */
                        gsl_matrix_set(*pmdPlm,il,2, -0.25*(il+2)*(il+1)*il*(il-1));
                }
            }
        }

        goto out;
    }

    xs = x*x-1;
    xr = sqrt(xs);

    /* DLMF 14.7.15; http://dlmf.nist.gov/14.7#ii */
    for(im = 1; im <= m; im++)
        gsl_matrix_set(mPlm, im, im, (2*im-1)*xr*gsl_matrix_get(mPlm, im-1, im-1));

    /* DLMF 14.10.7; http://dlmf.nist.gov/14.10 */
    for(im = 0; im <= MIN(m,l-1); im++)
        gsl_matrix_set(mPlm, im+1, im, (2*im+1)*x*gsl_matrix_get(mPlm,im,im));

    /* DLMF 14.10.3; http://dlmf.nist.gov/14.10 */
    for(im = 0; im <= m; im++)
        for(il = im+2; il <= l; il++)
            gsl_matrix_set(mPlm,il,im, ((2*il-1)*x* gsl_matrix_get(mPlm,il-1,im)-(im+il-1)*gsl_matrix_get(mPlm,il-2,im))/(il-im));

    /* calculate Plm' - if needed */
    if(pmdPlm != NULL)
    {
        /* DLMF 14.10.5; http://dlmf.nist.gov/14.10 */
        for(im = 1; im <= l; im++)
            gsl_matrix_set(*pmdPlm,im,0, im*(x*gsl_matrix_get(mPlm,im,0)-gsl_matrix_get(mPlm,im-1,0))/xs);

        /* derivative of DLMF 14.7.11 & DLMF 14.10.6 */
        for(im = 1; im <= m; im++)
            for(il=im; il <= l; il++)
                gsl_matrix_set(*pmdPlm,il,im, -im*x*gsl_matrix_get(mPlm,il,im)/xs+(il+im)*(il-im+1)/xr*gsl_matrix_get(mPlm,il,im-1));
    }

    out:
        result = gsl_matrix_get(mPlm, l, m);

        if(pmPlm != NULL)
            *pmPlm = mPlm;
        else
            gsl_matrix_free(mPlm);

        return result;
}

double dPlm(int l, int m, double x)
{
    double result;
    gsl_matrix *dplm;

    Plm(l,m,x,NULL,&dplm);
    result = gsl_matrix_get(dplm,l,m);
    gsl_matrix_free(dplm);

    return result;
}

double Yl1mYl2m(int l1, int l2, int m, double x)
{
    int lmax = MAX(l1,l2);
    gsl_matrix *plm = NULL;
    double result;
    
    Plm(lmax, m, x, &plm, NULL);
    result = -pow(-1, l2+m+(m%2))*Nlm(l1,m)*gsl_matrix_get(plm,l1,m)/sqrt(l1*(l1+1)) * Nlm(l2,m)*gsl_matrix_get(plm,l2,m)/sqrt(l2*(l2+1));
    gsl_matrix_free(plm);

    return result;
}

double Pl1mPl2m(int l1, int l2, int m, double x)
{
    int lmax = MAX(l1,l2);
    gsl_matrix *plm = NULL;
    double result;
    
    Plm(lmax, m, x, &plm, NULL);
    result = pow(-1, l2+m+(m%2))*gsl_matrix_get(plm,l1,m)*gsl_matrix_get(plm,l2,m);
    gsl_matrix_free(plm);

    return result;
}

double dYl1mdYl2m(int l1, int l2, int m, double x)
{
    int lmax = MAX(l1,l2);
    gsl_matrix *dplm = NULL;
    double result;
    
    Plm(lmax, m, x, NULL, &dplm);
    result = -pow(-1, l2+m+1-(m%2))*Nlm(l1,m)*gsl_matrix_get(dplm,l1,m)/sqrt(l1*(l1+1)) * Nlm(l2,m)*gsl_matrix_get(dplm,l2,m)/sqrt(l2*(l2+1));
    gsl_matrix_free(dplm);

    return result;
}

double dPl1mdPl2m(int l1, int l2, int m, double x)
{
    int lmax = MAX(l1,l2);
    gsl_matrix *dplm = NULL;
    double result;
    
    Plm(lmax, m, x, NULL, &dplm);
    result = pow(-1, l2+m+1-(m%2))*gsl_matrix_get(dplm,l1,m)*gsl_matrix_get(dplm,l2,m);
    gsl_matrix_free(dplm);

    return result;
}

double Yl1mdYl2m(int l1, int l2, int m, double x)
{
    int lmax = MAX(l1,l2);
    gsl_matrix *plm, *dplm = NULL;
    double result;
    
    Plm(lmax, m, x, &plm, &dplm);
    result = -pow(-1, l2+m+1-(m%2))*Nlm(l1,m)*gsl_matrix_get(plm,l1,m)/sqrt(l1*(l1+1)) * Nlm(l2,m)*gsl_matrix_get(dplm,l2,m)/sqrt(l2*(l2+1));
    gsl_matrix_free(plm);
    gsl_matrix_free(dplm);

    return result;
}

double Pl1mdPl2m(int l1, int l2, int m, double x)
{
    int lmax = MAX(l1,l2);
    gsl_matrix *plm, *dplm = NULL;
    double result;
    
    Plm(lmax, m, x, &plm, &dplm);
    result = pow(-1, l2+m+1-(m%2))*gsl_matrix_get(plm,l1,m)*gsl_matrix_get(dplm,l2,m);
    gsl_matrix_free(plm);
    gsl_matrix_free(dplm);

    return result;
}
