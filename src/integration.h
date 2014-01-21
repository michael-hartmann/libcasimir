#ifndef __INTEGRATION_H
#define __INTEGRATION_H

#include <quadmath.h>
#include "libcasimir.h"

void inline polymult(__float128 p1[], size_t len_p1, __float128 p2[], size_t len_p2, __float128 pdest[]);
double polyintegrate(__float128 p[], size_t len, int l1, int l2, int m, __float128 scale);
void polym(__float128 p[], int m, __float128 xi);
void polyplm(__float128 pl1[], __float128 pl2[], int l1, int l2, int m, __float128 xi);
void polydplm(__float128 pl1[], __float128 pl2[], int l1, int l2, int m, __float128 xi);
double log_polyintegrate(__float128 p[], size_t len, int l1, int l2, int m, int *sign);
void casimir_integrate(casimir_integrals_t *cint, int l1, int l2, int m, __float128 xi);

#endif
