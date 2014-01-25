#ifndef __INTEGRATION_H
#define __INTEGRATION_H

#include "libcasimir.h"

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
void casimir_integrate(casimir_integrals_t *cint, int l1, int l2, int m, double xi);

#endif
