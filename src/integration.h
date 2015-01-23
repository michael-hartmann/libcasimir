#ifndef __INTEGRATION_H
#define __INTEGRATION_H

#include "libcasimir.h"
#include "edouble.h"

typedef struct {
    int sign_A;
    edouble lnA_TE;
    edouble lnA_TM;

    int sign_B;
    edouble lnB_TE;
    edouble lnB_TM;
    
    int sign_C;
    edouble lnC_TE;
    edouble lnC_TM;
    
    int sign_D;
    edouble lnD_TE;
    edouble lnD_TM;
} integrands_drude_t;

void casimir_integrate_perf(casimir_integrals_t *cint, int l1, int l2, int m, double nT);
void casimir_integrate_drude(casimir_t *self, casimir_integrals_t *cint, int l1, int l2, int m, double nT);

void integrands_drude(edouble x, integrands_drude_t *integrands, casimir_t *self, double nT, int l1, int l2, int m);

void inline polymult(edouble p1[], size_t len_p1, edouble p2[], size_t len_p2, edouble pdest[]);
double log_polyintegrate(edouble p[], size_t len, int l1, int l2, int m, int *sign);
void polym(edouble p[], int m, edouble xi);
void polyplm(edouble pl1[], edouble pl2[], int l1, int l2, int m, edouble xi);
void polydplm(edouble pl1[], edouble pl2[], int l1, int l2, int m, edouble xi);

#endif
