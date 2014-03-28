#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libcasimir.h"

int main(int argc, char *argv[])
{
    casimir_t casimir;
    casimir_mie_cache_t cache;
    int n = 1;
    int m = 1;
    double Q = 0.99;
    double T = 0.1;
    int lmax = 100;
    double logdet = 0;

    casimir_init(&casimir, Q, T);
    casimir_set_lmax(&casimir, lmax);

    casimir_mie_cache_init(&cache, n*casimir.T*casimir.RbyScriptL);
    casimir_mie_cache_alloc(&casimir, &cache, casimir.lmax);

    logdet = casimir_logdetD(&casimir, n, m, &cache);
    casimir_free(&casimir);

    printf("logdet=%g\n", logdet);

    return 0;
}
