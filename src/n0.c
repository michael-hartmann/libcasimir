#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "casimir.h"
#include "libcasimir.h"

/* This function returns the seconds since 1st Jan 1970 in µs precision */
double now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec + tv.tv_usec*1e-6;
}

int main(int argc, char *argv[])
{
    double lfac = 6;
    int lmax, m, verbose = 1;
    double precision = 1e-8;
    double start_time = now();
    double Q, RbyL, F = 0;
    casimir_t casimir;

    if(argc < 2)
    {
        fprintf(stderr, "Usage: %s Q [lfac, precision]\n", argv[0]);
        return 1;
    }

    if(argc > 2)
        lfac = atof(argv[2]);
    if(argc > 3)
        precision = atof(argv[3]);

    Q = atof(argv[1]);
    RbyL = Q/(1-Q);
    lmax = ceil(RbyL*lfac);

    fprintf(stderr, "# Q=%g, R/L=%g, precision=%g, lfac=%g, lmax=%d\n", Q, RbyL, precision, lfac, lmax);
    
    casimir_init(&casimir, Q, 0.1);
    casimir_set_precision(&casimir, precision);
    casimir_set_verbose(&casimir, verbose);
    casimir_set_lmax(&casimir, lmax);

    for(m = 0; m < lmax; m++)
    {
        double value = casimir_logdetD(&casimir,0,m,NULL);

        value /= 2; // n == 0
        if(m == 0)
            value /= 2;

        fprintf(stderr, "# m=%d, value=%g\n", m, value);
        F += value;
        if(fabs(value/F) < precision)
            break;
    }
    
    casimir_free(&casimir);
    
    fprintf(stderr, "# Q, F, time\n");
    printf("%.15g, %.15g, %g\n", Q, F, now()-start_time);

    return 0;
}
