#ifndef __INTEGRATION_H
#define __INTEGRATION_H

    #include "libcasimir.h"

    void casimir_integrate(casimir_integrals_t *cint, int l1, int l2, int m, double nT);

        #include "edouble.h"

        void inline polymult(edouble p1[], size_t len_p1, edouble p2[], size_t len_p2, edouble pdest[]);
        double log_polyintegrate(edouble p[], size_t len, int l1, int l2, int m, int *sign);
        void polym(edouble p[], int m, edouble xi);
        void polyplm(edouble pl1[], edouble pl2[], int l1, int l2, int m, edouble xi);
        void polydplm(edouble pl1[], edouble pl2[], int l1, int l2, int m, edouble xi);
#endif
