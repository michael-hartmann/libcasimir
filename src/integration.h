#ifndef __INTEGRATION_H
#define __INTEGRATION_H

#include "libcasimir.h"

void inline polymult(double p1[], size_t len_p1, double p2[], size_t len_p2, double pdest[]);
double polyintegrate(double p[], size_t len, int l1, int l2, int m, double scale);
void polym(double p[], int m, double xi);
void polyplm(double pl1[], double pl2[], int l1, int l2, int m, double xi);
void polydplm(double pl1[], double pl2[], int l1, int l2, int m, double xi);
double log_polyintegrate(double p[], size_t len, int l1, int l2, int m, int *sign);
void casimir_integrate(casimir_integrals_t *cint, int l1, int l2, int m, double xi);

#endif
