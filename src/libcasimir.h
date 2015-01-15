#ifndef __LIBCASIMIR_H
#define __LIBCASIMIR_H

#include <pthread.h>

#define HBARC 3.161526510740123e-26
#define KB    1.3806488e-23

#define CASIMIR_DEFAULT_PRECISION 1e-12
#define CASIMIR_IDLE 1000
#define CASIMIR_FACTOR_LMAX 5

const char *casimir_compile_info(void);

/**
 * The Casimir object. This structure stores all essential information about
 * temperature, geometry and the reflection properties of the mirrors.
 *
 * Do not modify the attributes of the structure yourself!
 */
typedef struct
{
    /**
     * @name geometry and temperature
     */
     /*@{*/
    double RbyScriptL; /**< \f$R/\mathcal{L}\f$, where \f$R\f$ is the radius of the sphere and \f$L\f$ is the separation of plane and sphere. */
    double T; /**< temperature */
    /*@}*/

    /**
     * @name reflection properties of the mirrors
     */
     /*@{*/
    double omegap_sphere; /**< plasma frequency \f$\omega_\mathrm{P}\f$ of sphere */
    double omegap_plane;  /**< plasma frequency \f$\omega_\mathrm{P}\f$ of plane */
    double gamma_sphere;  /**< relaxation frequency \f$\gamma\f$ of sphere */
    double gamma_plane;   /**< relaxation frequency \f$\gamma\f$ of plane */
    /*@}*/

    /**
     * @name accuracy and numerical parameters
     */
     /*@{*/
    int lmax;            /**< truncation value for vector space \f$\ell_\mathrm{max}\f$ */
    int verbose;         /**< flag that indicates to be verbose */
    int extrapolate;     /**< flag that indicates to use extrapolation */
    int cores;           /**< number of thread that should be used */
    double precision;    /**< precision */
    pthread_t **threads; /**< list of pthread objects */
    /*@}*/
} casimir_t;


/**
 * thread object.
 */
typedef struct
{
    casimir_t *self; /**< pointer to Casimir object */
    int n;           /**< Matsubara term */
    int nmax;        /**< maximum number of n */
    double value;    /**< free energy for Matsubara term n*/
} casimir_thread_t;



/**
 * Cache for Mie coefficients.
 */
typedef struct
{
    double *al;   /**< list of Mie coefficients \f$a_\ell\f$ (logarithms) */
    int *al_sign; /**< list of signs of Mie coefficients \f$a_\ell\f$ */
    double *bl;   /**< list of Mie coefficients \f$b_\ell\f$ (logarithms) */
    int *bl_sign; /**< list of signs of Mie coefficients \f$b_\ell\f$ */
    int lmax;     /**< truncation value for vector space \f$\ell_\mathrm{max}\f$ */
    int n;        /**< Matsubara term */
} casimir_mie_cache_t;


/*
typedef struct
{
    double logA,logB,logC,logD;
    int signA, signB, signC, signD;
} casimir_integrals_t;
*/

typedef struct
{
    double lnA_TE, lnA_TM;
    double lnB_TE, lnB_TM;
    double lnC_TE, lnC_TM;
    double lnD_TE, lnD_TM;
    int signA_TE, signA_TM;
    int signB_TE, signB_TM;
    int signC_TE, signC_TM;
    int signD_TE, signD_TM;
} casimir_integrals_t;


/* prototypes */
double casimir_epsilon(double xi, double omegap, double gamma_);
double casimir_lnepsilon(double xi, double omegap, double gamma_);

double casimir_lnLambda(int l1, int l2, int m, int *sign);
double casimir_lnXi(int l1, int l2, int m, int *sign);

double casimir_F_SI_to_scaled(double F_SI, double ScriptL_SI);
double casimir_F_scaled_to_SI(double F, double ScriptL_SI);
double casimir_T_SI_to_scaled(double T_SI, double ScriptL_SI);
double casimir_T_scaled_to_SI(double T, double ScriptL_SI);

void casimir_lnab0(int l, double *a0, int *sign_a0, double *b0, int *sign_b0);

int casimir_init(casimir_t *self, double RbyScriptL, double T);
void casimir_free(casimir_t *self);

int casimir_set_omegap_sphere(casimir_t *self, double omegap);
double casimir_get_omegap_sphere(casimir_t *self);
int casimir_set_gamma_sphere(casimir_t *self, double gamma_);
double casimir_get_gamma_sphere(casimir_t *self);


int casimir_get_lmax(casimir_t *self);
int casimir_set_lmax(casimir_t *self, int lmax);

int casimir_get_cores(casimir_t *self);
int casimir_set_cores(casimir_t *self, int cores);

double casimir_get_precision(casimir_t *self);
int    casimir_set_precision(casimir_t *self, double precision);

int casimir_get_verbose(casimir_t *self);
int casimir_set_verbose(casimir_t *self, int verbose);

int casimir_get_extrapolate(casimir_t *self);
int casimir_set_extrapolate(casimir_t *self, int extrapolate);

void casimir_lnab(casimir_t *self, const int n, const int l, double *lna, double *lnb, int *sign_a, int *sign_b);
double casimir_lna_perf(casimir_t *self, const int l, const int n, int *sign);
double casimir_lnb_perf(casimir_t *self, const int l, const int n, int *sign);

double casimir_F_n(casimir_t *self, const int n, int *mmax);
double casimir_F(casimir_t *self, int *nmax);

void casimir_mie_cache_init(casimir_mie_cache_t *cache, int n);
int casimir_mie_cache_alloc(casimir_t *self, casimir_mie_cache_t *cache);
void casimir_mie_cache_free(casimir_mie_cache_t *cache);

double casimir_logdetD0(casimir_t *self, int m, double *EE, double *MM);
double casimir_logdetD(casimir_t *self, int n, int m, casimir_mie_cache_t *cache);

void casimir_rp(casimir_t *self, double nT, double k, double *r_TE, double *r_TM);

#endif
