#define _ISOC99_SOURCE
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

#include "la.h"
#include "logdet1m.h"
#include "plm_fast.h"
#include "libcasimir.h"
#include "gausslaguerre.h"

#define FACTOR_LMAX 3

/*
 * This is the error handler for the Gnu scientifiy library.
 * In case there are any errors, print a debug message and exit
 */
void handler(const char *reason, const char *file, int line, int gsl_errno)
{
    fprintf(stderr, "# %s in %s:%d (gsl_errno: %d)\n", reason, file, line, gsl_errno);
    exit(1);
}

/* Lambda 
 * This function returns the Λ prefactor for given l1,l2,m.
 *
 * The values are computed using the lngamma function, so that nominator
 * and denominator of the term
 * (l1-m)!*(l2-m)!/(l1+m)!/(l2+m) 
 * don't overflow.
 *
 * Restrictions: l1,l2,m integers, l1,l2>=1, l1,l2 >= m
 * Symmetries: Λ(l1,l2,m) = Λ(l2,l1,m)
 */
double Lambda(int l1, int l2, int m)
{
    return -sqrt( ((double)2*l1+1)*(2*l2+1)/(l1*l2*(l1+1)*(l2+1)) ) \
         * sqrt( exp(gsl_sf_lngamma(l1-m+1)+gsl_sf_lngamma(l2-m+1)-gsl_sf_lngamma(l1+m+1)-gsl_sf_lngamma(l2+m+1)) );
}


/* casimir_logdet1m
 * This function returns returns log(det(1-M))
 *
 * The value is calculated either using Mercator series, if |M| < 1 or
 * by calculating the eigenvalues of M
 */
int casimir_logdet1m(gsl_matrix *M, double *logdet, int n, int m, const char *desc)
{
    int ret;
    double norm = la_norm_froebenius(M);

    // if |M| < 1 use Mercator series, otherwise use eigenvalues
    /*
    if(norm <= 0.97)
        ret = logdet1m_taylor(M, logdet);
    else
        ret = logdet1m_eigenvalues(M, logdet);
    */

    la_matrix_info(stderr, M);
    ret = logdet1m_eigenvalues(M, logdet);
    if(ret != 0)
    {
        int i,j;
        int nans = 0, infs = 0, zeros = 0;
        for(i = 0; i < M->size1; i++)
            for(j = 0; j < M->size1; j++)
            {
                double elem = gsl_matrix_get(M, i,j);
                if(isinf(elem)) infs++;
                if(isnan(elem)) nans++;
                if(elem == 0)   zeros++;
            }

        fprintf(stderr, "# Can't calculate log(det(1-%s)): %d (n=%d, m=%d, norm=%g, %d zeros, %d nans, %d infs, dim %dx%d)\n", \
                        desc, ret, n, m, norm, zeros, nans, infs, (int)M->size1, (int)M->size2);
    }

    return ret;
}

/* Function Xi
 * This function returns the Ξ prefactor for given l1,l2,m
 *
 * lgamma is used to prevent overflows - like in Λ.
 *
 * Restrictions: l1,l2,m integers, l1,l2>=1, l1,l2 >= m
 */
double Xi(int l1,int l2, int m)
{
    return pow(-1, l2)/pow(4, 2*l1+l2+1) \
         * exp(log(-Lambda(l1,l2,m))+gsl_sf_lngamma(2*l1+1)+gsl_sf_lngamma(2*l2+1)+gsl_sf_lngamma(l1+l2+1)-gsl_sf_lngamma(l1)-gsl_sf_lngamma(l2)-gsl_sf_lngamma(l1-m+1)-gsl_sf_lngamma(l2-m+1) );
}

/*
 * Initializes casimir object for perfect mirrors.
 * See casimir_init
 */
int casimir_init_perfect(casimir_t *self, double RbyL, double T)
{
    return casimir_init(self, RbyL, T, INFINITY, 0);
}

/*
 * The Casimir class provides methods to calculate the free energy of the
 * Casimir effect in the plane-sphere geometry for perfect reflectors.
 *
 * The main goal is to calculate the free energy and derived quantities.
 * Create a new Casimir object.
 * 
 * R/L:    ratio of radius of sphere and distance between sphere and plate
 * T:      Temperature
 * omegap: plasma frequency for Drude of plasma model
 * gamma:  relaxation frequency for Drude modell
 *
 * Restrictions: RbyL,T,omegap,gamma > 0
 */
int casimir_init(casimir_t *self, double RbyL, double T, double omegap, double gamma)
{
    if(RbyL < 0 || T < 0)
        return -1;

    gsl_set_error_handler(&handler);

    /*
    self->int_limits    = INT_LIMITS;
    self->int_workspace = gsl_integration_workspace_alloc(self->int_limits);
    self->epsrel        = INT_EPSREL;
    */

    self->lmax = (int)ceil(RbyL*FACTOR_LMAX);

    self->T       = T;
    self->RbyL    = RbyL;
    self->omegap  = omegap;
    self->gamma   = gamma;
    self->eps_n   = 1e-9;
    self->verbose = 0;

    return 0;
}

/*
 * Set maximum value for l
 */
void casimir_set_lmax(casimir_t *self, int lmax)
{
    self->lmax = lmax;
}

void casimir_set_limits(casimir_t *self, int limits)
{
    self->int_limits = limits;
}

void casimir_set_epsrel(casimir_t *self, double epsrel)
{
    self->epsrel = epsrel;
}

void casimir_set_verbose(casimir_t *self, int verbose)
{
    self->verbose = verbose;
}

void casimir_set_eps_n(casimir_t *self, double eps_n)
{
    self->eps_n = eps_n;
}

/*
 * Free memory for casimir object
 */
void casimir_free(casimir_t *self)
{
    //gsl_integration_workspace_free(self->int_workspace);
}


/*
 * For small x<<1 a_l will scale as
 * 
 * a_l(x) ~ a0*(x/2)^(2l+1)
 *
 * This method returns the prefactor a0
 */
double a0(int l)
{
    return M_PI*pow(-1, l)*( 2*gsl_sf_gamma(l+1.5)-l*gsl_sf_gamma(l+0.5) )/( l*pow(gsl_sf_gamma(l+0.5),2)*gsl_sf_gamma(l+1.5) );
}

/*
 * For small x<<1 b_l will scale as
 * 
 * b_l(x) ~ b0*(x/2)^(2l+1)
 *
 * This method returns the prefactor b0
 */
double b0(int l)
{
    return M_PI*pow(-1, l+1)/( gsl_sf_gamma(l+0.5)*gsl_sf_gamma(l+1.5) );
}

/*
 * Returns the coefficient a_l for reflection on the sphere
 *
 * Restrictions: l integer, l>=1, xi>0
 */
double casimir_a(casimir_t *self, int l, double arg)
{
    //double n,n2;

    if(isinf(self->omegap))
        return M_PI/2*pow(-1,l+1)*(l*gsl_sf_bessel_Inu(l+0.5,arg)-arg*gsl_sf_bessel_Inu(l-0.5,arg))/(l*gsl_sf_bessel_Knu(l+0.5,arg)+arg*gsl_sf_bessel_Knu(l-0.5,arg));

    /*
    n2 = epsilon(self,arg);
    n  = sqrt(n2);

    return M_PI/2*pow(-1, l+1)*(n2*SLA(l,n,arg)-SLB(l,n,arg))/(n2*SLC(l,n,arg)-SLD(l,n,arg));
    */

    return 0;
}

/*        
 * Returns the coefficient b_l for reflection on the sphere
 *
 * Restrictions: l integer, l>=1, xi>0
 */        
double casimir_b(casimir_t *self, int  l, double arg)
{
    //double n;

    if(isinf(self->omegap))
        /* perfect reflector */
        return M_PI/2*pow(-1, l+1)*gsl_sf_bessel_Inu(l+0.5,arg)/gsl_sf_bessel_Knu(l+0.5,arg);

    /*
    n = sqrt(epsilon(self,arg));

    return M_PI/2*pow(-1, l+1)*(SLA(l,n,arg)-SLB(l,n,arg))/(SLC(l,n,arg)-SLD(l,n,arg));
    */
    return 0;
}

/*
 * Returns dielectric function epsilon(omega)
 */
 /*
double epsilon(casimir_t *self, double xi)
{
    return 1+gsl_pow_2(self->omegap)/(xi*(xi+self->gamma));
}
*/

/*
 * Returns the TE reflection coefficient for the plate
 */
double r_TE(casimir_t *self, double x, double arg)
{
    //double cos2, root, eps;

    if(isinf(self->omegap))
        /* perfect reflector */
        return -1;
    
    /*
    cos2 = gsl_pow_2(1 + x*self->c/(2*self->L*xi)); // = kappa*c/xi
    eps  = epsilon(self, xi);
    root = sqrt(1+(eps-1)/cos2);

    return (1-root)/(1+root);
    */

    return 0;
}

/*
 * Returns the TM reflection coefficient for the plate
 */
double r_TM(casimir_t *self, double x, double xi)
{
    //double cos2, root, eps;

    if(isinf(self->omegap))
        /* perfect reflector */
        return +1;

    /*
    cos2 = gsl_pow_2(1 + x*self->c/(2*self->L*xi));
    eps  = epsilon(self,xi);
    root = sqrt(1+(eps-1)/cos2);

    return (eps-root)/(eps+root);
    */

    return 0;
}

/*
 * Returns the integrands A,B,C,D for l1,l2,m at the point x. len should be
 * always 4. The values of A(x), B(x), C(x) and D(x) are stored in vec:
 *     vec[] = { A(x), B(x), C(x), D(x) }
 * This function is called by gausslaguerre_integrate_vec.
 */
void casimir_integrands_vec(double x, void *params, double vec[4], int len)
{
    casimir_int_t *p = (casimir_int_t *)params;
    int l1 = p->l1;
    int l2 = p->l2;
    int m  = p->m;
    double nT = p->nT;
    // fac = -2 / sqrt( l1*(l1+1) * l2*(l2+1) ) * r_p
    double fac = -2/sqrt(l1*l2*(l1+1)*(l2+1)) * p->rp(p->self,x,nT) * p->scale;
    double Yl1m, Yl2m, dYl1m, dYl2m;

    /* Note, that we split up Lambda here:
     *     Lambda = -2/sqrt( l1*(l1+1) * l2*(l2+1) ) * Nl1m * Nl2m
     * Nl1m and Nl2m is in [d]Yl[1,2]m, the factor 
     *     -2/sqrt( l1*(l1+1) * l2*(l2+1) ) * Nl1m * Nl2m
     * is in fac.
     *
     * scale is just a factor, that may be used to get finite results.
     * 
     * To integrate f(x)*exp(x), we just need f(x). Therefore we don't need the
     * factor exp(x) here.
     */

    plm_Yl12md(l1, l2, m, 1+x/(2*nT), &Yl1m, &Yl2m, &dYl1m, &dYl2m);

    vec[0] = fac * pow(-1, l2+m+(m%2))   * Yl1m  * Yl2m  * 1/(gsl_pow_2(x)+4*nT*x);
    vec[1] = fac * pow(-1, l2+m+1-(m%2)) * dYl1m * dYl2m *   (gsl_pow_2(x)+4*nT*x);
    vec[2] = fac * pow(-1, l2+m-(m%2))   * Yl1m  * dYl2m;
    vec[3] = fac * pow(-1, l1+m-(m%2))   * Yl2m  * dYl1m;
}

/*
 * Returns the integrals A,B,C,D for l1,l2,m,xi and p=TE,TM
 */
int casimir_integrate(casimir_t *self, casimir_integrals_t *cint, int l1, int l2, int n, int m, double scale)
{
    double vec[4];
    double prefactor[4];
    double nT = n*self->T;
    int N = ceil((l1+l2-m+3)/2.);
    casimir_int_t params;

    prefactor[0] = exp(-2*nT)*gsl_pow_2(m)*2*nT;
    prefactor[1] = exp(-2*nT)/gsl_pow_3(2*nT);
    prefactor[2] = exp(-2*nT)*m/(2*nT);
    prefactor[3] = pow(-1, l1+l2+1)*prefactor[2];

    params.self  = self;
    params.l1    = l1;
    params.l2    = l2;
    params.m     = m;
    params.nT    = nT;
    params.rp    = r_TM;
    params.scale = scale;

    /* integrate */
    gausslaguerre_integrate_vec(casimir_integrands_vec, &params, N, vec, 4);
    cint->A_TM = prefactor[0]*vec[0];
    cint->B_TM = prefactor[1]*vec[1];
    cint->C_TM = prefactor[2]*vec[2];
    cint->D_TM = prefactor[3]*vec[3];

    /* if we have perfect mirrors, we don't need to calculate the
     * TE integrals as they just differ in sign */
    if(isinf(self->omegap))
    {
        cint->A_TE = -cint->A_TM;
        cint->B_TE = -cint->B_TM;
        cint->C_TE = -cint->C_TM;
        cint->D_TE = -cint->D_TM;
    }
    else
    {
        params.rp = r_TE;
        gausslaguerre_integrate_vec(casimir_integrands_vec, &params, n, vec, 4);
        cint->A_TE = prefactor[0]*vec[0];
        cint->B_TE = prefactor[1]*vec[1];
        cint->C_TE = prefactor[2]*vec[2];
        cint->D_TE = prefactor[3]*vec[3];
    }

    return 0;
}


/*
 * Allocate memory for the Mie-coefficients a_l and b_l
 */
int casimir_mie_cache_alloc(casimir_t *self, casimir_mie_cache_t *cache, double arg)
{
    int l1, lmax = self->lmax;

    if(arg == 0)
    {
        cache->al = cache->bl = NULL;
        return 0;
    }

    cache->al = (double *)malloc((lmax+1)*sizeof(double));
    cache->bl = (double *)malloc((lmax+1)*sizeof(double));

    if(cache->al == NULL || cache->bl == NULL)
    {
        fprintf(stderr, "# Out of memory.\n");
        exit(1);
    }

    cache->al[0] = cache->bl[0] = 0;
    for(l1 = 1; l1 <= lmax; l1++)
    {
        cache->al[l1] = casimir_a(self,l1,arg);
        cache->bl[l1] = casimir_b(self,l1,arg);
    }

    return 1;
}

/*
 * Free memory of cache.
 */
void casimir_mie_cache_free(casimir_mie_cache_t *cache)
{
    if(cache->al != NULL)
        free(cache->al);
    if(cache->bl != NULL)
        free(cache->bl);

    cache->al = cache->bl = NULL;
}

/*
 * Calculate free energy. Sum from nmin to nmax.

 * Restrictions: nmax integer, nmax >= 0
 */
double casimir_F(casimir_t *self, int *nmax)
{
    int n = 0;
    double sum0, sum = 0;
    double RbyL = self->RbyL;

    while(1)
    {
        casimir_mie_cache_t cache;
        double sum_n = 0;
        int m;

        casimir_mie_cache_alloc(self, &cache, n*self->T*RbyL);

        for(m = 0; m <= self->lmax; m++)
        {
            double value = casimir_logdetD(self,n,m,&cache);

            if(m == 0)
                value /= 2;

            sum_n += value;
        }

        if(n == 0)
        {
            sum0 = sum_n;
            sum += sum_n/2;
        }
        else
            sum += sum_n;

        casimir_mie_cache_free(&cache);

        //fprintf(stderr, "# %d, %.15g\n", n, sum_n);
        if(sum_n/sum0 < self->eps_n)
        {
            if(nmax != NULL)
                *nmax = n;
            return self->T/M_PI*sum;
        }

        n++;
    }
}

/*
 * Calculate logarithm of determinant of D=1-M for given m,ξ
 *
 * Restrictions: m integer, m>=0, ξ>= 0
 */
double casimir_logdetD(casimir_t *self, int n, int m, casimir_mie_cache_t *cache)
{
    int min,max,dim,l1,l2;
    double logdet_EE = 0;
    double logdet_MM = 0;
    double logdet = 0;

    min = MAX(m,1);
    max = self->lmax;
    
    dim = (max-min+1);
    
    if(n == 0)
    {
        gsl_matrix *EE, *MM;
        double RbyL = self->RbyL;

        EE = gsl_matrix_alloc(dim, dim);
        MM = gsl_matrix_alloc(dim, dim);

        for(l1 = min; l1 <= max; l1++)
            for(l2 = min; l2 <= max; l2++)
            {
                double XiRL = Xi(l1,l2,m)*pow(RbyL, l1+l2+1);

                gsl_matrix_set(EE, l1-min, l2-min, +a0(l1)*XiRL);
                gsl_matrix_set(MM, l1-min, l2-min, -b0(l1)*XiRL);
            }
    
        casimir_logdet1m(EE, &logdet_EE, n, m, "M_EE");
        casimir_logdet1m(MM, &logdet_MM, n, m, "M_MM");

        gsl_matrix_free(EE);
        gsl_matrix_free(MM);

        return logdet_EE+logdet_MM;
    }
    else
    {
        gsl_matrix *M = gsl_matrix_alloc(2*dim, 2*dim);
    
        /* M_EE, -M_EM
           M_ME,  M_MM */
        for(l1 = min; l1 <= max; l1++)
        {
            for(l2 = min; l2 <= l1; l2++)
            {
                casimir_integrals_t cint;
                double scale = 1;
                double nTRbyL = n*self->T*self->RbyL;
                double al1 = cache->al[l1];
                double bl1 = cache->bl[l1];
                double al2 = cache->al[l2];
                double bl2 = cache->bl[l2];

                if(nTRbyL < 1)
                {
                    double scale_l1, scale_l2;
                    scale    = pow(nTRbyL, +l2+l1);
                    scale_l2 = pow(nTRbyL, 2*l1);
                    scale_l1 = pow(nTRbyL, 2*l2);

                    al1 /= scale_l1;
                    al2 /= scale_l2;
                    bl1 /= scale_l1;
                    bl2 /= scale_l2;

                    if(al1 == 0)
                        al1 = a0(l1)/pow(2, l1+1)*pow(nTRbyL, 1-l1);
                    if(al2 == 0)
                        al2 = a0(l2)/pow(2, l2+1)*pow(nTRbyL, 1-l2);
                    if(bl1 == 0)
                        bl1 = b0(l1)/pow(2, l1+1)*pow(nTRbyL, 1-l1);
                    if(bl2 == 0)
                        bl2 = b0(l2)/pow(2, l2+1)*pow(nTRbyL, 1-l2);
                }

                casimir_integrate(self, &cint, l1, l2, n, m, scale);

                gsl_matrix_set(M,     l1-min,     l2-min, al1*(cint.A_TE+cint.B_TM)); /* M_EE */
                gsl_matrix_set(M,     l2-min,     l1-min, pow(-1, l1+l2)*al2*(cint.A_TE+cint.B_TM)); /* M_EE */

                gsl_matrix_set(M, dim+l1-min, dim+l2-min, bl1*(cint.A_TM+cint.B_TE)); /* M_MM */
                gsl_matrix_set(M, dim+l2-min, dim+l1-min, pow(-1, l1+l2)*bl2*(cint.A_TM+cint.B_TE)); /* M_MM */

                gsl_matrix_set(M, dim+l1-min,     l2-min, al1*(cint.C_TE+cint.D_TM)); /* M_EM */
                gsl_matrix_set(M, dim+l2-min,     l1-min, pow(-1, l1+l2+1)*al2*(cint.D_TE+cint.C_TM)); /* M_EM */

                gsl_matrix_set(M,     l1-min, dim+l2-min, bl1*(cint.C_TM+cint.D_TE)); /* - M_ME */
                gsl_matrix_set(M,     l2-min, dim+l1-min, pow(-1, l1+l2+1)*bl2*(cint.D_TM+cint.C_TE)); /* - M_ME */
            }
        }

        if(m == 0)
        {
            int i,j;
            double logdet_EE = 0, logdet_MM = 0;
            gsl_matrix *EE = gsl_matrix_alloc(dim, dim);
            gsl_matrix *MM = gsl_matrix_alloc(dim, dim);

            for(i = 0; i < dim; i++)
                for(j = 0; j < dim; j++)
                {
                    gsl_matrix_set(EE, i,j, gsl_matrix_get(M, i,j));
                    gsl_matrix_set(MM, i,j, gsl_matrix_get(M, dim+i,dim+j));
                }

            casimir_logdet1m(EE, &logdet_EE, n, m, "EE");
            casimir_logdet1m(MM, &logdet_MM, n, m, "MM");

            logdet = logdet_EE+logdet_MM;

            gsl_matrix_free(EE);
            gsl_matrix_free(MM);
        }
        else
            casimir_logdet1m(M, &logdet, n, m, "M");

        gsl_matrix_free(M);
        return logdet;
    }
}
