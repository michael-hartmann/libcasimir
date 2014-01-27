#ifndef __INTEGRATION_H
#define __INTEGRATION_H

    #include "libcasimir.h"

    void casimir_integrate(casimir_integrals_t *cint, int l1, int l2, int m, double xi);

    #ifdef INTEGRATION_QUAD
        #include "quad.h"

        void inline polymult(__float128 p1[], size_t len_p1, __float128 p2[], size_t len_p2, __float128 pdest[]);
        double log_polyintegrate(__float128 p[], size_t len, int l1, int l2, int m, int *sign);
        void polym(__float128 p[], int m, __float128 xi);
        void polyplm(__float128 pl1[], __float128 pl2[], int l1, int l2, int m, __float128 xi);
        void polydplm(__float128 pl1[], __float128 pl2[], int l1, int l2, int m, __float128 xi);
    #else
        typedef struct {
            double value;
            int sign;
        } log_t;

        void polyprint(log_t p[], size_t len);

        void inline polymult(log_t p1[], size_t len_p1, log_t p2[], size_t len_p2, log_t pdest[]);
        double polyintegrate(log_t p[], size_t len, int l1, int l2, int m, double scale);
        void polym(log_t p[], int m, double xi);
        void polyplm(log_t pl1[], log_t pl2[], int l1, int l2, int m, double xi);
        void polydplm(log_t pl1[], log_t pl2[], int l1, int l2, int m, double xi);
        double log_polyintegrate(log_t p[], size_t len, int l1, int l2, int m, int *sign);
    #endif

#endif
