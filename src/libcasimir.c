/**
 * @file   libcasimir.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   October, 2014
 * @brief  library to calculate the free Casimir energy in the plane-sphere geometry
 */


#define _GNU_SOURCE

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <pthread.h>
#include <unistd.h>

#include "edouble.h"
#include "integration.h"
#include "libcasimir.h"
#include "matrix.h"
#include "sfunc.h"
#include "utils.h"


static char CASIMIR_COMPILE_INFO[255];

/** @brief Return string with compile information
 *
 * The returned string contains information which compiler was used and which
 * kind of arithmetics the binary uses.
 *
 * Do not modify this string!
 *
 * @retval constant string
 */
const char *casimir_compile_info(void)
{
    snprintf(CASIMIR_COMPILE_INFO, sizeof(CASIMIR_COMPILE_INFO)/sizeof(char), "Compiler %s, using %s", COMPILER, CASIMIR_ARITHMETICS);
    return CASIMIR_COMPILE_INFO;
}

void casimir_info(casimir_t *self, FILE *stream, const char *prefix)
{
    if(prefix == NULL)
        prefix = "";

    fprintf(stream, "%sRbyScriptL = %g\n", prefix, self->RbyScriptL);
    fprintf(stream, "%sT = %g\n", prefix, self->T);

    fprintf(stream, "%somegap_sphere   = %g\n", prefix, self->omegap_sphere);
    fprintf(stream, "%somegap_plane    = %g\n", prefix, self->omegap_plane);
    fprintf(stream, "%sgamma_sphere    = %g\n", prefix, self->gamma_sphere);
    fprintf(stream, "%sgamma_plane     = %g\n", prefix, self->gamma_plane);
    fprintf(stream, "%sintegration     = ", prefix);
    if(self->integration <= 0)
        fprintf(stream, "analytic (perfect mirrors)\n");
    else
        fprintf(stream, "%d\n", self->integration);

    fprintf(stream, "%slmax = %d\n",        prefix,  self->lmax);
    fprintf(stream, "%sverbose = %d\n",     prefix,  self->verbose);
    fprintf(stream, "%sextrapolate = %d\n", prefix, self->extrapolate);
    fprintf(stream, "%scores = %d\n",       prefix, self->cores);
    fprintf(stream, "%sprecision = %g\n",   prefix, self->precision);
}


/**
 * @brief Calculate logarithm and sign of prefactor \f$\Lambda_{\ell_1 \ell_2}^{(m)}\f$
 *
 * This function returns the logarithm of the prefactor for given
 * \f$\ell_1,\ell_2,m\f$. This prefactor is defined by (cf Eq. (5.19))
 * \f[
 *      \Lambda_{\ell_1,\ell_2}^{(m)} = -\frac{2 N_{\ell_1,m} N_{\ell_2,m}}{\sqrt{\ell_1 (\ell_1+1) \ell_2 (\ell_2+1)}}
 * \f]
 *
 * If sign is not NULL, -1 will be stored in sign.
 *
 * The values are computed using the lngamma function in a smart way to avoid overflows.
 *
 * Restrictions: \f$\ell_1,\ell_2 \ge 1\f$, \f$\ell_1,\ell_2 \ge m\f$
 *
 * Symmetries: \f$\Lambda_{\ell_1,\ell_2}^{(m)} = \Lambda_{\ell_2,\ell_1}^{(m)}\f$
 *
 * @param [in]  l1
 * @param [in]  l2
 * @param [in]  m
 * @param [out] sign
 * @retval log(Lambda(l1,l2,m))
 */
double inline casimir_lnLambda(int l1, int l2, int m, int *sign)
{
    if(sign != NULL)
        *sign = -1;
    return LOG2 + (logq(2.*l1+1)+logq(2*l2+1)-LOG4-logq(l1)-logq(l1+1)-logq(l2)-logq(l2+1)+lnfac(l1-m)+lnfac(l2-m)-lnfac(l1+m)-lnfac(l2+m))/2.0L;
}


/**
 * @brief Calculate \f$\epsilon(i\xi)\f$ for Drude model
 *
 * This function returns the dielectric function
 * \f[
 *      \epsilon(i\xi) = 1 + \frac{\omega_\mathrm{P}^2}{\xi(\xi+\gamma)}
 * \f]
 *
 * @param [in]  xi     imaginary frequency (in scaled units: \f$\xi=nT\f$)
 * @param [in]  omegap Plasma frequency
 * @param [in]  gamma_ relaxation frequency
 * @retval epsilon(xi, omegap, gamma_)
 */
double casimir_epsilon(double xi, double omegap, double gamma_)
{
    return 1+ omegap*omegap/(xi*(xi+gamma_));
}


/**
 * @brief Calculate \f$\log \epsilon(i\xi)\f$ for Drude model
 *
 * This function returns the logarithm of the dielectric function
 * \f[
 *      \epsilon(i\xi) = 1 + \frac{\omega_\mathrm{P}^2}{\xi(\xi+\gamma)}
 * \f]
 *
 * @param [in]  xi     imaginary frequency (in scaled units: \f$\xi=nT\f$)
 * @param [in]  omegap Plasma frequency
 * @param [in]  gamma_ relaxation frequency
 * @retval log(epsilon(xi, omegap, gamma_))
 */
double casimir_lnepsilon(double xi, double omegap, double gamma_)
{
    return log1p(omegap*omegap/(xi*(xi+gamma_)));
}


/**
 * @brief Calculate Fresnel coefficients \f$r_{TE}\f$ and \f$r_{TM}\f$ for Drude model
 *
 * This function calculates the Fresnel coefficients for TE and TM mode
 *
 * @param [in]      self    Casimir object
 * @param [in]      nT      imaginary frequency (in scaled units: \f$\xi=nT\f$)
 * @param [in]      k       xy projection of wavevector
 * @param [in,out]  r_TE    Fresnel coefficient for TE mode
 * @param [in,out]  r_TM    Fresnel coefficient for TM mode
 */
void casimir_rp(casimir_t *self, double nT, double k, double *r_TE, double *r_TM)
{
    double epsilon = casimir_epsilon(nT, self->omegap_plane, self->gamma_plane);
    double beta = sqrt(1 + (epsilon-1)/(1 + pow_2(k/nT)));

    *r_TE = (1-beta)/(1+beta);
    *r_TM = (epsilon-beta)/(epsilon+beta);
}

/**
* @name converting
*/
/*@{*/

/**
 * @brief Convert free energy in SI units to free energy in units of \f$\mathcal{L}/(\hbar c)\f$
 *
 * This function returns 
 * \f[
 *      \mathcal{F}_\mathrm{scaled} = \mathcal{F}_\mathrm{SI} \frac{\mathcal{L}}{\hbar c}
 * \f]
 *
 * @param [in] F_SI free energy in SI units
 * @param [in] ScriptL \f$\mathcal{L} = R+L\f$ (in units of meters)
 * @retval free energy in scaled units
 */
double casimir_F_SI_to_scaled(double F_SI, double ScriptL)
{
    return ScriptL/(HBARC)*F_SI;
}


/**
 * @brief Convert free energy in units of \f$\mathcal{L}/(\hbar c)\f$ to free energy in SI units
 *
 * This function returns 
 * \f[
 *      \mathcal{F}_\mathrm{SI} = \mathcal{F}_\mathrm{scaled} \frac{\hbar c}{\mathcal{L}}
 * \f]
 *
 * @param [in] F free energy in units of \f$\mathcal{L}/(\hbar c)\f$
 * @param [in] ScriptL \f$\mathcal{L} = R+L\f$ (in units of meters)
 * @retval free energy in SI units
 */
double casimir_F_scaled_to_SI(double F, double ScriptL)
{
    return HBARC/ScriptL*F;
}


/**
 * @brief Convert temperature in units of Kelvin to temperature in units of \f$2\pi k_B \mathcal{L}/(\hbar c)\f$
 *
 * This function returns 
 * \f[
 *      T_\mathrm{scaled} = \frac{2\pi k_b \mathcal{L}}{\hbar c} T_\mathrm{SI}
 * \f]
 *
 * @param [in] T_SI temperature in units of Kelvin
 * @param [in] ScriptL \f$\mathcal{L} = R+L\f$ (in units of meters)
 * @retval temperature in unitss of \f$2\pi k_B \mathcal{L}/(\hbar c)\f$
 */
double casimir_T_SI_to_scaled(double T_SI, double ScriptL)
{
    return 2*PI*KB*ScriptL/HBARC*T_SI;
}


/**
 * @brief Convert temperature in units of \f$2\pi k_B \mathcal{L}/(\hbar c)\f$ to temperature in units of Kelvin
 *
 * This function returns 
 * \f[
 *      T_\mathrm{scaled} = \frac{2\pi k_b \mathcal{L}}{\hbar c} T_\mathrm{SI}
 * \f]
 *
 * @param [in] T temperature in units of \f$2\pi k_B \mathcal{L}/(\hbar c)\f$
 * @param [in] ScriptL \f$\mathcal{L} = R+L\f$ (in units of meters)
 * @retval temperature in units of Kelvin
 */
double casimir_T_scaled_to_SI(double T, double ScriptL)
{
    return HBARC/(2*PI*KB*ScriptL)*T;
}

/*@}*/

/**
 * @brief Calculate logarithm and sign of prefactor \f$\Xi_{\ell_1 \ell_2}^{(m)}\f$
 *
 * This function returns the logarithm of the prefactor for given
 * \f$\ell_1,\ell_2,m\f$. The prefactor is defined by Eq. (5.54).
 *
 * If sign is not NULL, the sign of \f$\Xi_{\ell_1 \ell_2}^{(m)}\f$ is stored in
 * sign.
 *
 * The values are computed using the lngamma function in a smart way to avoid overflows.
 *
 * Restrictions: \f$\ell_1,\ell_2 \ge 1\f$, \f$\ell_1,\ell_2 \ge m\f$
 *
 * @param [in]  l1
 * @param [in]  l2
 * @param [in]  m
 * @param [out] sign
 * @retval log(Xi(l1,l2,m))
 */
double casimir_lnXi(int l1, int l2, int m, int *sign)
{
    if(sign != NULL)
        *sign = pow(-1, l2);
    return (log(2*l1+1)+log(2*l2+1)-lnfac(l1-m)-lnfac(l2-m)-lnfac(l1+m)-lnfac(l2+m)-log(l1)-log(l1+1)-log(l2)-log(l2+1))/2.0 \
           +lnfac(2*l1)+lnfac(2*l2)+lnfac(l1+l2)-LOG4*(2*l1+l2+1)-lnfac(l1-1)-lnfac(l2-1);
}

/**
* @name initialization and setting parameters
*/
/*@{*/

/**
 * @brief Create a new Casimir object for perfect reflectors
 *
 * Restrictions: \f$T > 0\f$, \f$0 < \mathcal{L} < 1\f$
 *
 * @param [out] self Casimir object
 * @param [in]  RbyScriptL \f$\frac{R}{\mathcal{L}} = \frac{R}{R+L}\f$
 * @param [in]  T temperature in units of \f$2\pi k_B \mathcal{L}/(\hbar c)\f$
 */
int casimir_init(casimir_t *self, double RbyScriptL, double T)
{
    double LbyR = 1./RbyScriptL - 1;
    if(RbyScriptL < 0 || RbyScriptL >= 1 || T < 0)
        return 0;

    self->lmax = (int)ceil(CASIMIR_FACTOR_LMAX/LbyR);

    self->T           = T;
    self->RbyScriptL  = RbyScriptL;
    self->precision   = CASIMIR_DEFAULT_PRECISION;
    self->extrapolate = 0;
    self->verbose     = 0;
    self->cores       = 1;
    self->threads     = NULL;
    
    /* perfect reflectors */
    self->integration = -1; /* perfect reflectors */
    self->omegap_sphere = INFINITY;
    self->gamma_sphere  = 0;
    self->omegap_plane  = INFINITY;
    self->gamma_plane   = 0;

    return 1;
}


/**
 * @brief Set order of integration
 *
 * Set order/type of integration.
 *
 * @param [in,out] self Casimir object
 * @param [in] integration: 0 perfect reflectors, >0: order of Gauss-Laguerre integration
 */
void casimir_set_integration(casimir_t *self, int integration)
{
    if(integration <= 0)
        self->integration = 0;
    else
        self->integration = integration;
}

/**
 * @brief Get order of integration
 *
 * Get order/type of integration.
 *
 * @param [in,out] self Casimir object
 * @retval 0 if analytic integration for perfect reflectors
 * @retval order of integration for Gauss-Laguerre
 */
int casimir_get_integration(casimir_t *self)
{
    return self->integration;
}


/**
 * @brief Set \f$\omega_\mathrm{P}\f$ for the sphere
 *
 * Set the plasma frequency for the sphere.
 *
 * @param [in,out] self Casimir object
 * @param [in] omegap plasma frequency
 * @retval 1 if successful
 * @retval 0 if omegap < 0
 */
int casimir_set_omegap_sphere(casimir_t *self, double omegap)
{
    if(omegap > 0)
    {
        self->omegap_sphere = omegap;
        self->integration   = 50;
        return 1;
    }
    return 0;
}

/**
 * @brief Set \f$\omega_\mathrm{P}\f$ for the plane
 *
 * Set the plasma frequency for the plane.
 *
 * @param [in,out] self Casimir object
 * @param [in] omegap plasma frequency
 * @retval 1 if successful
 * @retval 0 if omegap < 0
 */
int casimir_set_omegap_plane(casimir_t *self, double omegap)
{
    if(omegap > 0)
    {
        self->omegap_plane = omegap;
        self->integration  = 50;
        return 1;
    }
    return 0;
}


/**
 * @brief Get \f$\omega_\mathrm{P}\f$ for the sphere
 *
 * Get the plasma frequency for the sphere.
 *
 * @param [in,out] self Casimir object
 * @retval plasma frequency
 */
double casimir_get_omegap_sphere(casimir_t *self)
{
    return self->omegap_sphere;
}

/**
 * @brief Get \f$\omega_\mathrm{P}\f$ for the plane
 *
 * Get the plasma frequency for the plane.
 *
 * @param [in,out] self Casimir object
 * @retval plasma frequency
 */
double casimir_get_omegap_plane(casimir_t *self)
{
    return self->omegap_plane;
}


/**
 * @brief Set \f$\gamma\f$ for the sphere
 *
 * Set the relaxation frequency for the sphere.
 *
 * @param [in,out] self Casimir object
 * @param [in] gamma_ relaxation frequency
 * @retval 1 if successful
 * @retval 0 if gamma_ < 0
 */
int casimir_set_gamma_sphere(casimir_t *self, double gamma_)
{
    if(gamma_ > 0)
    {
        self->gamma_sphere = gamma_;
        self->integration = 50;
        return 1;
    }
    return 0;
}

/**
 * @brief Set \f$\gamma\f$ for the plane
 *
 * Set the relaxation frequency for the plane.
 *
 * @param [in,out] self Casimir object
 * @param [in] gamma_ relaxation frequency
 * @retval 1 if successful
 * @retval 0 if gamma_ < 0
 */
int casimir_set_gamma_plane(casimir_t *self, double gamma_)
{
    if(gamma_ > 0)
    {
        self->gamma_plane = gamma_;
        self->integration = 50;
        return 1;
    }
    return 0;
}


/**
 * @brief Get \f$\gamma\f$ for the sphere
 *
 * Get the relaxation frequency for the sphere.
 *
 * @param [in,out] self Casimir object
 * @retval relaxation frequency
 */
double casimir_get_gamma_sphere(casimir_t *self)
{
    return self->gamma_sphere;
}

/**
 * @brief Get \f$\gamma\f$ for the plane
 *
 * Get the relaxation frequency for the plane.
 *
 * @param [in,out] self Casimir object
 * @retval relaxation frequency
 */
double casimir_get_gamma_plane(casimir_t *self)
{
    return self->gamma_plane;
}


/**
 * @brief Get extrapolation flag
 *
 * Extrapolation is experimental and considered dangerous at the moment.
 *
 * @param [in,out] self Casimir object
 * @retval Extrapolation flag
 */
int casimir_get_extrapolate(casimir_t *self)
{
    return self->extrapolate;
}


/**
 * @brief Set extrapolation flag
 *
 * Extrapolation is experimental and considered dangerous at the moment.
 *
 * @param [in,out] self Casimir object
 * @param extrapolate extrapolation flag
 * @retval 1
 */
int casimir_set_extrapolate(casimir_t *self, int extrapolate)
{
    self->extrapolate = extrapolate ? 1 : 0;
    return 1;
}


/**
 * @brief Return numbers of used cores
 *
 * See casimir_set_cores.
 *
 * @param [in,out] self Casimir object
 * @retval number of used cores (>=0)
 */
int casimir_get_cores(casimir_t *self)
{
    return self->cores;
}


/**
 * @brief Set the number of used cores
 *
 * This library supports multiple processor cores. However, you must specify
 * how many cores the library should use. By default, only one core will be
 * used. If you have a quad core computer, you might want to set the number of
 * cores to 4.
 *
 * The libraray uses POSIX threads for parallelization. Threads share memory and
 * for this reason all cores must be on the same computer.
 *
 * Restrictions: cores > 0
 *
 * @param [in,out] self Casimir object
 * @param [in] cores number of cores that should be used
 * @retval 1 if successful
 * @retval 0 if cores < 1
 */
int casimir_set_cores(casimir_t *self, int cores)
{
    if(cores < 1)
        return 0;

    self->cores = cores;
    self->threads = xrealloc(self->threads, cores*sizeof(pthread_t));

    return 1;
}


/**
 * @brief Set maximum value of l
 *
 * In general the round trip matrices are infinite. For a numerical evaluation
 * the dimension has to be limited to a finite value. The accuracy of the
 * result depends on the truncation of the vector space. For more information,
 * cf. chapter 6.1.
 * 
 * @param [in,out] self Casimir object
 * @param [in] lmax maximum number of l
 * @retval 1 if successful
 * @retval 0 if lmax < 1
 */
int casimir_set_lmax(casimir_t *self, int lmax)
{
    if(lmax <= 0)
        return 0;

    self->lmax = lmax;
    return 1;
}


/**
 * @brief Get maximum value of l
 *
 * See casimir_set_lmax.
 *
 * @param [in,out] self Casimir object
 * @retval lmax maximum value of l
 */
int casimir_get_lmax(casimir_t *self)
{
    return self->lmax;
}


/**
 * @brief Get verbose flag
 *
 * Return if the verbose flag is set.
 *
 * @param [in,out] self Casimir object
 * @retval 0 if verbose flag is not set
 * @retval 1 if verbose flag is set
 */
int casimir_get_verbose(casimir_t *self)
{
    return self->verbose;
}


/**
 * @brief Set verbose flag
 *
 * Use this function to set the verbose flag. If set to 1, this will cause the
 * library to print information to stderr. 
 *
 * @param [in,out] self Casimir object
 * @param [in] verbose 1 if verbose, 0 if not verbose
 * @retval 1
 */
int casimir_set_verbose(casimir_t *self, int verbose)
{
    self->verbose = verbose ? 1 : 0;
    return 1;
}


/**
 * @brief Get precision
 *
 * See casimir_set_precision
 *
 * @param [in,out] self Casimir object
 * @retval precision
 */
double casimir_get_precision(casimir_t *self)
{
    return self->precision;
}


/**
 * @brief Set precision
 *
 * @param [in,out] self Casimir object
 * @param [in] precision
 * @retval 1 if successful
 * @retval 0 if precision <= 0
 */
int casimir_set_precision(casimir_t *self, double precision)
{
    if(precision <= 0)
        return 0;

    self->precision = precision;
    return 1;
}


/*
 * Free memory for casimir object
 */
/**
 * @brief Free memory for Casimir object
 *
 * This function will free allocated memory for the Casimir object. If you have
 * allocated memory for the object yourself, you have, however, to free this
 * yourself.
 * 
 * @param [in,out] self Casimir object
 */
void casimir_free(casimir_t *self)
{
    if(self->threads != NULL)
    {
        xfree(self->threads);
        self->threads = NULL;
    }
}

/*@}*/


/**
* @name Mie coefficients
*/
/*@{*/

/** Return the logarithm of the prefactors \f$a_{\ell,0}^\mathrm{perf}\f$, \f$b_{\ell,0}^\mathrm{perf}\f$ and its signs
 *
 * For small frequencies \f$\chi = \frac{\xi R}{c} \ll 1\f$ the Mie
 * coeffiecients scale like
 * \f[
 * a_{\ell}^\mathrm{perf} = a_{\ell,0}^\mathrm{perf} \left(\frac{\chi}{2}\right)^{2\ell+1} \\
 * \f]
 * \f[
 * b_{\ell}^\mathrm{perf} = b_{\ell,0}^\mathrm{perf} \left(\frac{\chi}{2}\right)^{2\ell+1}
 * \f]
 * This function returns the logarithm of the prefactors
 * \f$a_{\ell,0}^\mathrm{perf}\f$, \f$b_{\ell,0}^\mathrm{perf}\f$ and its
 * signs.
 *
 * In scaled units: \f$\chi = nT \frac{R}{\mathcal{L}}\f$
 *
 * @param [in] l
 * @param [out] a0 coefficient \f$a_{\ell,0}^\mathrm{perf}\f$
 * @param [out] sign_a0 sign of \f$a_{\ell,0}^\mathrm{perf}\f$
 * @param [out] b0 coefficient \f$b_{\ell,0}^\mathrm{perf}\f$
 * @param [out] sign_b0 sign of \f$b_{\ell,0}^\mathrm{perf}\f$
 */
void casimir_lnab0(int l, double *a0, int *sign_a0, double *b0, int *sign_b0)
{
    *sign_a0 = pow(-1, l);
    *sign_b0 = pow(-1, l+1);
    *b0 = LOGPI-lngamma(l+0.5)-lngamma(l+1.5);
    *a0 = *b0+log1p(1.0/l);
}


/**
 * @brief Return logarithm of Mie coefficient \f$a_\ell\f$ for perfect reflectors and its sign
 *
 * The frequency will be determined by n: \f$\xi = nT\f$
 *
 * Restrictions: \f$\ell \ge 1\f$, \f$\ell \ge 0\f$
 *
 * @param [in,out] self Casimir object
 * @param [in] l
 * @param [in] n Matsubara term, \f$xi = nT\f$
 * @param [out] sign sign of \f$a_\ell\f$
 * @retval logarithm of Mie coefficient a_l
 */
double casimir_lna_perf(casimir_t *self, const int l, const int n, int *sign)
{
    edouble nominator, denominator, frac, ret;
    edouble lnKlp,lnKlm,lnIlm,lnIlp;
    edouble prefactor;
    edouble chi = n*self->T*self->RbyScriptL;
    edouble lnfrac = log(chi)-log(l);

    /* we could do both calculations together. but it doesn't cost much time -
     * so why bother? 
     */
    bessel_lnInuKnu(l-1, chi, &lnIlm, &lnKlm);
    bessel_lnInuKnu(l,   chi, &lnIlp, &lnKlp);

    prefactor = LOGPI-LOG2+lnIlp-lnKlp;
    *sign = pow(-1, l+1);

    /* numinator */
    {
        frac = expq(lnfrac+lnIlm-lnIlp);
        if(frac < 1)
            nominator = log1pq(fabsq(-frac));
        else
        {
            if(frac > 1)
                *sign *= -1;

            nominator = logq(fabsq(1-frac));
        }
    }
    /* denominator */
    {
        frac = expq(lnfrac+lnKlm-lnKlp);
        if(frac < 1)
            denominator = log1pq(frac);
        else
            denominator = log1pq(frac);
    }

    ret = prefactor+nominator-denominator;

    assert(!isnan(ret));
    assert(!isinf(ret));

    return ret;
}


/**
 * @brief Return logarithm of Mie coefficient \f$b_\ell\f$ for perfect reflectors and its sign
 *
 * The frequency will be determined by n: \f$\xi = nT\f$
 *
 * Restrictions: \f$\ell \ge 1\f$, \f$\ell \ge 0\f$
 *
 * @param [in,out] self Casimir object
 * @param [in] l
 * @param [in] n Matsubara term, \f$xi = nT\f$
 * @param [out] sign sign of \f$b_\ell\f$
 * @retval logarithm of Mie coefficient b_l
 */
double casimir_lnb_perf(casimir_t *self, const int l, const int n, int *sign)
{
    edouble chi = n*self->T*self->RbyScriptL;
    edouble lnInu, lnKnu, ret;

    bessel_lnInuKnu(l, chi, &lnInu, &lnKnu);
    *sign = pow(-1, l+1);

    ret = LOGPI-LOG2+lnInu-lnKnu;

    assert(!isnan(ret));
    assert(!isinf(ret));

    return ret;
}


/**
 * @brief Return logarithm of Mie coefficients \f$a_\ell\f$, \f$b_\ell\f$ for Drude model
 *
 * For \f$\omega_\mathrm{P} = \infty\f$ the Mie coefficient for perfect
 * reflectors are returned (see casimir_lna_perf and casimir_lnb_perf).
 *
 * Cf. Eqs. (3.30) and (3.31).
 *
 * sign_a and sign_b must be valid pointers and must not be NULL.
 *
 * @param [in,out] self Casimir object
 * @param [in] n Matsubara term, \f$\xi = nT\f$
 * @param [in] l
 * @param [out] lna logarithm of Mie coefficient \f$a_\ell\f$
 * @param [out] lnb logarithm of Mie coefficient \f$b_\ell\f$
 * @param [out] sign_a sign of Mie coefficient \f$a_\ell\f$
 * @param [out] sign_b sign of Mie coefficient \f$b_\ell\f$
 */
void casimir_lnab(casimir_t *self, const int n_mat, const int l, double *lna, double *lnb, int *sign_a, int *sign_b)
{ 
    int sign_sla, sign_slb, sign_slc, sign_sld;
    edouble ln_n, ln_sla, ln_slb, ln_slc, ln_sld;
    edouble lnIl, lnKl, lnIlm, lnKlm, lnIl_nchi, lnKl_nchi, lnIlm_nchi, lnKlm_nchi;
    edouble xi = n_mat*self->T;
    edouble chi = xi*self->RbyScriptL;
    edouble ln_chi = logq(xi)+logq(self->RbyScriptL);
    edouble omegap = self->omegap_sphere;
    edouble gamma_ = self->gamma_sphere;
    int sign_a_num, sign_a_denom, sign_b_num, sign_b_denom;

    if(isinf(omegap))
    {
        *lna = casimir_lna_perf(self, l, n_mat, sign_a);
        *lnb = casimir_lnb_perf(self, l, n_mat, sign_b);
        return;
    }

    ln_n = casimir_lnepsilon(xi, omegap, gamma_)/2;

    bessel_lnInuKnu(l,   chi, &lnIl,  &lnKl);
    bessel_lnInuKnu(l-1, chi, &lnIlm, &lnKlm);

    bessel_lnInuKnu(l,   expq(ln_n)*chi, &lnIl_nchi,  &lnKl_nchi);
    bessel_lnInuKnu(l-1, expq(ln_n)*chi, &lnIlm_nchi, &lnKlm_nchi);

    ln_sla = lnIl_nchi + logadd_s(lnIl,      +1, ln_chi+lnIlm,           -1, &sign_sla);
    ln_slb = lnIl      + logadd_s(lnIl_nchi, +1, ln_n+ln_chi+lnIlm_nchi, -1, &sign_slb);
    ln_slc = lnIl_nchi + logadd_s(lnKl,      +1, ln_chi+lnKlm,           +1, &sign_slc);
    ln_sld = lnKl      + logadd_s(lnIl_nchi, +1, ln_n+ln_chi+lnIlm_nchi, -1, &sign_sld);

    /*
    printf("n2=%.14g\n", n2);
    printf("lnIl = %.14g\n", exp(lnIl));
    printf("\n");
    printf("chi=%.14g\n", chi);
    */

    /*
    printf("sla=%.14g\n", sign_sla*exp(ln_sla));
    printf("slb=%.14g\n", sign_slb*exp(ln_slb));
    printf("slc=%.14g\n", sign_slc*exp(ln_slc));
    printf("sld=%.14g\n", sign_sld*exp(ln_sld));
    */

    /* XXX calculate this somehow in a smarter way */
    *lna = LOGPI - LOG2 + logadd_s(2*ln_n+ln_sla, +sign_sla, ln_slb, -sign_slb, &sign_a_num) - logadd_s(2*ln_n+ln_slc, +sign_slc, ln_sld, -sign_sld, &sign_a_denom);
    *lnb = LOGPI - LOG2 + logadd_s(       ln_sla, +sign_sla, ln_slb, -sign_slb, &sign_b_num) - logadd_s(       ln_slc, +sign_slc, ln_sld, -sign_sld, &sign_b_denom);

    *sign_a = sign_a_num*sign_a_denom;
    *sign_b = sign_b_num*sign_b_denom;

    assert(!isnan(*lna));
    assert(!isinf(*lna));
    assert(!isnan(*lnb));
    assert(!isinf(*lnb));
}


/**
 * @brief Initialize Mie cache
 *
 * The cache stores all Mie coefficients \f$a_\ell\f$ and \f$b_\ell\f$ for the
 * imaginary frequency \f$\xi = nT\f$. This function initialized the cache.
 * However, no Mie coefficients will be computed and stored.
 *
 * @param [out] cache cache for Mie coefficients
 * @param [in] n Matsubara term, \f$\xi=nT\f$
 */
void casimir_mie_cache_init(casimir_mie_cache_t *cache, int n)
{
    cache->al = cache->bl = NULL;
    cache->al_sign = cache->bl_sign = NULL;
    cache->lmax = 0;
    cache->n = n;
}


/**
 * @brief Allocate memory for the Mie coefficients \f$a_\ell\f$ and \f$b_\ell\f$
 *
 * This function computes all the Mie coefficients for \f$\xi=nT\f$ and stores
 * it in cache. Make sure you have already initialized the cache (cf.
 * casimir_mie_cache_init).
 *
 * @param [in,out] self Casimir object
 * @param [in,out] cache Mie cache
 */
int casimir_mie_cache_alloc(casimir_t *self, casimir_mie_cache_t *cache)
{
    int n = cache->n;
    int l, lmax = self->lmax;

    if(n == 0)
    {
        cache->al = cache->bl = NULL;
        cache->al_sign = cache->bl_sign = NULL;
        return 0;
    }

    cache->al      = (double *)xrealloc(cache->al,      (lmax+1)*sizeof(double));
    cache->bl      = (double *)xrealloc(cache->bl,      (lmax+1)*sizeof(double));
    cache->al_sign =    (int *)xrealloc(cache->al_sign, (lmax+1)*sizeof(int));
    cache->bl_sign =    (int *)xrealloc(cache->bl_sign, (lmax+1)*sizeof(int));

    cache->al[0] = cache->bl[0] = 0;
    for(l = MAX(1,cache->lmax); l <= lmax; l++)
        casimir_lnab(self, n, l, &cache->al[l], &cache->bl[l], &cache->al_sign[l], &cache->bl_sign[l]);

    cache->lmax = lmax;

    return 1;
}

/**
 * @brief Free memory of cache.
 *
 * This function will free allocated memory for the cache.
 *
 * @param [in, out] cache Mie cache
 */
void casimir_mie_cache_free(casimir_mie_cache_t *cache)
{
    if(cache->al != NULL)
        xfree(cache->al);
    if(cache->bl != NULL)
        xfree(cache->bl);
    if(cache->al_sign != NULL)
        xfree(cache->al_sign);
    if(cache->bl_sign != NULL)
        xfree(cache->bl_sign);

    cache->al = cache->bl = NULL;
}

/*@}*/

/* Sum len numbers in value.
   The idea is: To avoid a loss of significance, we sum beginning with smallest
   number and add up in increasing order
*/
static double _sum(double values[], size_t len)
{
    int i;
    double sum = 0;

    for(i = len-1; i > 0; i--)
        sum += values[i];

    sum += values[0]/2;

    return sum;
}

static void *_thread_f(void *p)
{
    casimir_thread_t *r = (casimir_thread_t *)p;
    r->value = casimir_F_n(r->self, r->n, &r->nmax);
    return r;
}

static pthread_t *_start_thread(casimir_thread_t *r)
{
    pthread_t *t = xmalloc(sizeof(pthread_t));
    pthread_create(t, NULL, _thread_f, (void *)r);

    return t;
}

static int _join_threads(casimir_t *self, double values[], int *ncalc)
{
    int i, joined = 0, running = 0;
    casimir_thread_t *r;
    pthread_t **threads = self->threads;

    for(i = 0; i < self->cores; i++)
    {
        if(threads[i] != NULL)
        {
            running++;

            if(pthread_tryjoin_np(*threads[i], (void *)&r) == 0)
            {
                joined++;

                if(r->n > *ncalc)
                    *ncalc = r->n;

                values[r->n] = r->value;
                xfree(r);
                xfree(threads[i]);
                threads[i] = NULL;

                if(self->verbose)
                    fprintf(stderr, "# n=%d, value=%.15g\n", r->n, values[r->n]);
            }
        }
    }

    if(running == 0)
        return -1;

    return joined;
}


/**
* @name Calculate free energy
*/
/*@{*/

/**
 * @brief Calculate free energy for Matsubara term n
 *
 * This function calculates the free energy for the Matsubara term n. If mmax
 * is not NULL, the maximum usedd value of m is stored in mmax.
 *
 * @param [in,out] self Casimir object
 * @param [in] n Matsubara term, \f$\xi=nT\f$
 * @param [out] mmax maximum number of m
 * @retval Casimir free energy for given n
 */
double casimir_F_n(casimir_t *self, const int n, int *mmax)
{
    double precision = self->precision;
    casimir_mie_cache_t cache;
    double sum_n = 0;
    int m;
    const int lmax = self->lmax;
    double values[lmax+1];

    for(m = 0; m <= lmax; m++)
        values[m] = 0;

    casimir_mie_cache_init(&cache, n);
    casimir_mie_cache_alloc(self, &cache);

    for(m = 0; m <= self->lmax; m++)
    {
        values[m] = casimir_logdetD(self,n,m,&cache);

        if(self->verbose)
            fprintf(stderr, "# n=%d, m=%d, value=%.15g\n", n, m, values[m]);

        /* If F is !=0 and value/F < 1e-16, then F+value = F. The addition
         * has no effect.
         * As for larger m value will be even smaller, we can skip the
         * summation here. 
         */
        sum_n = _sum(values, lmax+1);
        if(values[0] != 0 && fabs(values[m]/sum_n) < precision)
            break;
    }

    casimir_mie_cache_free(&cache);

    if(self->verbose)
        fprintf(stderr, "# n=%d, value=%.15g\n", n, sum_n);

    if(mmax != NULL)
        *mmax = m;

    return sum_n;
}


/**
 * @brief Calculate free energy
 *
 * This function calculates the free energy. If nmax is not NULL, the highest
 * Matsubara term used will be stored in nnmax.
 * 
 * @param [in,out] self Casimir object
 * @param [out] nmax maximum number of n
 * @retval Casimir free energy
 */
double casimir_F(casimir_t *self, int *nmax)
{
    int i, n = 0;
    double sum_n = 0;
    const double precision = self->precision;
    double *values = NULL;
    size_t len = 0;
    int ncalc = 0;
    const int cores = self->cores;
    pthread_t **threads = self->threads;

    if(cores > 1)
        for(i = 0; i < cores; i++)
            threads[i] = NULL;

    /* So, here we sum up all m and n that contribute to F.
     * So, what do we do here?
     *
     * We want to evaluate
     *      \sum_{n=0}^\infty \sum{m=0}^{l_max} log(det(1-M)),
     * where the terms for n=0 and m=0 are weighted with a factor 1/2.
     */
    while(1)
    {
        if(n >= len)
        {
            const int delta = MAX(512, cores);

            values = (double *)xrealloc(values, (len+delta)*sizeof(double));

            for(i = len; i < len+delta; i++)
                values[i] = 0;

            len += delta;
        }

        if(cores > 1)
        {
            casimir_thread_t *r;

            for(i = 0; i < cores; i++)
            {
                if(threads[i] == NULL)
                {
                    r = (casimir_thread_t *)xmalloc(sizeof(casimir_thread_t));

                    r->self  = self;
                    r->n     = n++;
                    r->value = 0;
                    r->nmax  = 0;

                    threads[i] = _start_thread(r);
                }
            }

            if(_join_threads(self, values, &ncalc) == 0)
                usleep(CASIMIR_IDLE);
        }
        else
        {
            values[n] = casimir_F_n(self, n, NULL);

            ncalc = n;
            n++;
        }

        if(values[0] != 0)
        {
            if(fabs(values[ncalc]/(2*values[0])) < precision)
            {
                if(cores > 1)
                    while(_join_threads(self, values, &ncalc) != -1)
                        usleep(CASIMIR_IDLE);

                sum_n = _sum(values, n);
                if(self->extrapolate && n > 20)
                {
                    double r1 = values[n-1]/values[n-2];
                    double r2 = values[n-2]/values[n-3];
                    double r3 = values[n-3]/values[n-4];
                    double r4 = values[n-4]/values[n-5];
                    double r5 = values[n-5]/values[n-6];
                    double r  = (r1+r2+r3+r4+r5)/5;
                    sum_n += values[n-1]*r/(1-r);
                }
                /* get out of here */
                if(nmax != NULL)
                    *nmax = n-1; // we calculated n term from n=0,...,nmax=n-1

                if(values != NULL)
                    xfree(values);

                return self->T/PI*sum_n;
            }
        }
    }
}


/**
 * @brief Calculate \f$\log\det \mathcal{D}(\xi=0)\f$
 *
 * This function calculates the logarithm of the determinant of the scattering
 * operator D for the Matsubara term \f$n=0\f$.
 *
 * @param [in,out] self Casimir object
 * @param [in] m
 * @param [out] logdet_EE
 * @param [out] logdet_MM
 * @retval log det D(xi=0)
 */
double casimir_logdetD0(casimir_t *self, int m, double *logdet_EE, double *logdet_MM)
{
    int l1,l2,min,max,dim;
    double lnRbyScriptL = log(self->RbyScriptL);
    double value_EE, value_MM;

    min = MAX(m,1);
    max = self->lmax;

    dim = (max-min+1);

    matrix_edouble_t *EE = matrix_edouble_alloc(dim);
    matrix_edouble_t *MM = matrix_edouble_alloc(dim);

    /* calculate the logarithm of the matrix elements of D */
    for(l1 = min; l1 <= max; l1++)
        for(l2 = min; l2 <= max; l2++)
        {
            /* i: row of matrix, j: column of matrix */
            const int i = l1-min, j = l2-min;
            int sign_a0, sign_b0, sign_xi;
            double lna0, lnb0;
            double lnXiRL = casimir_lnXi(l1,l2,m,&sign_xi)+(2*l1+1)*lnRbyScriptL;
            casimir_lnab0(l1, &lna0, &sign_a0, &lnb0, &sign_b0);

            matrix_set(EE, i,j, (l1 == l2 ? 1 : 0) - sign_xi*sign_a0*expq(lna0+lnXiRL));
            matrix_set(MM, i,j, (l1 == l2 ? 1 : 0) - sign_xi*sign_a0*expq(lnb0+lnXiRL));
        }

    /* balance the matrix */
    matrix_edouble_balance(EE);
    matrix_edouble_balance(MM);

    value_EE = matrix_edouble_logdet(EE);
    value_MM = matrix_edouble_logdet(MM);

    /* free space for matrices */
    matrix_edouble_free(EE);
    matrix_edouble_free(MM);

    if(logdet_EE != NULL)
        *logdet_EE = value_EE;
    if(logdet_MM != NULL)
        *logdet_MM = value_MM;

    return value_EE+value_MM;
}


/**
 * @brief Calculate \f$\log\det \mathcal{D}(\xi=nT)\f$
 *
 * This function calculates the logarithm of the determinant of the scattering
 * operator D for the Matsubara term \f$n\f$.
 *
 * @param [in,out] self Casimir object
 * @param [in] n Matsubara term
 * @param [in] m
 * @param [in] cache Mie cache
 * @retval log det D(xi=nT)
 */
double casimir_logdetD(casimir_t *self, int n, int m, casimir_mie_cache_t *cache)
{
    int min,max,dim,l1,l2;
    double logdet = 0;
    double nTRbyScriptL = n*self->T*self->RbyScriptL;

    min = MAX(m,1);
    max = self->lmax;

    dim = (max-min+1);

    if(n == 0)
        return casimir_logdetD0(self, m, NULL, NULL);

    matrix_edouble_t *M = matrix_edouble_alloc(2*dim);

    /* M_EE, -M_EM
       M_ME,  M_MM */
    for(l1 = min; l1 <= max; l1++)
    {
        for(l2 = min; l2 <= l1; l2++)
        {
            int Delta_ij = (l1 == l2 ? 1 : 0);
            const int i = l1-min, j = l2-min;
            casimir_integrals_t cint;
            double lnal1 = cache->al[l1];
            double lnbl1 = cache->bl[l1];
            double lnal2 = cache->al[l2];
            double lnbl2 = cache->bl[l2];

            double al1_sign = cache->al_sign[l1];
            double bl1_sign = cache->bl_sign[l1];
            double al2_sign = cache->al_sign[l2];
            double bl2_sign = cache->bl_sign[l2];

            if(nTRbyScriptL < 1)
            {
                double lognTRbyScriptL = log(nTRbyScriptL);
                lnal1 -= (l1-l2)*lognTRbyScriptL;
                lnbl1 -= (l1-l2)*lognTRbyScriptL;

                lnal2 -= (l2-l1)*lognTRbyScriptL;
                lnbl2 -= (l2-l1)*lognTRbyScriptL;
            }

            if(self->integration > 0)
                casimir_integrate_drude(self, &cint, l1, l2, m, n*self->T);
            else
                casimir_integrate_perf(&cint, l1, l2, m, n*self->T);

            /* EE */
            matrix_set(M, i,j, Delta_ij -                al1_sign*( cint.signA_TE*expq(lnal1+cint.lnA_TE) + cint.signB_TM*expq(lnal1+cint.lnB_TM) ));
            matrix_set(M, j,i, Delta_ij - pow(-1, l1+l2)*al2_sign*( cint.signA_TE*expq(lnal2+cint.lnA_TE) + cint.signB_TM*expq(lnal2+cint.lnB_TM) ));

            /* MM */
            matrix_set(M, i+dim,j+dim, Delta_ij -                bl1_sign*( cint.signA_TM*expq(lnbl1+cint.lnA_TM) + cint.signB_TE*expq(lnbl1+cint.lnB_TE) ));
            matrix_set(M, j+dim,i+dim, Delta_ij - pow(-1, l1+l2)*bl2_sign*( cint.signA_TM*expq(lnbl2+cint.lnA_TM) + cint.signB_TE*expq(lnbl2+cint.lnB_TE) ));

            if(m != 0)
            {
                /* M_EM */
                matrix_set(M, dim+i,j, -                  al1_sign*( cint.signC_TE*expq(lnal1+cint.lnC_TE) + cint.signD_TM*expq(lnal1+cint.lnD_TM) ));
                matrix_set(M, dim+j,i, - pow(-1, l1+l2+1)*al2_sign*( cint.signD_TE*expq(lnal2+cint.lnD_TE) + cint.signC_TM*expq(lnal2+cint.lnC_TM) ));

                /* M_ME */
                matrix_set(M, i,dim+j, -                  bl1_sign*( cint.signC_TM*expq(lnbl1+cint.lnC_TM) + cint.signD_TE*expq(lnbl1+cint.lnD_TE) ));
                matrix_set(M, j,dim+i, - pow(-1, l1+l2+1)*bl2_sign*( cint.signD_TM*expq(lnbl2+cint.lnD_TM) + cint.signC_TE*expq(lnbl2+cint.lnC_TE) ));
            }
        }
    }

    if(m == 0)
    {
        size_t i,j;
        matrix_edouble_t *EE = matrix_edouble_alloc(dim);
        matrix_edouble_t *MM = matrix_edouble_alloc(dim);

        for(i = 0; i < dim; i++)
            for(j = 0; j < dim; j++)
            {
                matrix_set(EE, i,j, matrix_get(M, i,j));
                matrix_set(MM, i,j, matrix_get(M, dim+i,dim+j));
            }

        matrix_edouble_balance(MM);
        matrix_edouble_balance(EE);

        logdet = matrix_edouble_logdet(EE)+matrix_edouble_logdet(MM);

        matrix_edouble_free(EE);
        matrix_edouble_free(MM);
    }
    else
    {
        matrix_edouble_balance(M);
        logdet = matrix_edouble_logdet(M);
    }

    matrix_edouble_free(M);

    assert(!isinf(logdet));
    return logdet;
}

/*@}*/
