#include <ctype.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "sfunc.h"
#include "libcasimir.h"

#define PRECISION 1e-10

/* This function returns the seconds since 1st Jan 1970 in µs precision */
static double now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec + tv.tv_usec*1e-6;
}

/* print usage */
static void usage(FILE *stream)
{
    fprintf(stream, "Usage: casimir_logdetD [OPTIONS]\n\
This program will calculate the free Casimir energy for the plane-sphere \n\
geometry for given n,m,T,L/R. \n\
\n\
Mandatory options:\n\
    -Q Radius of sphere divided by distance between center of sphere and plate; 0 < R/(R+L) < 1\n\
    -T Temperature\n\
    -n value of n\n\
    -m value of m\n\
\n\
Further options:\n\
    -l, --lscale\n\
        Specify parameter lscale. The vector space has to be truncated at some\n\
        value lmax. This program will use lmax=(R/L*lscale) (default: 3)\n\
\n\
    -L\n\
        Set lmax to given value. When -L is used, -l will be ignored\n\
\n\
    -p, --precision\n\
        Set precision to given value. Default: %g\n\
\n\
    --buffering\n\
        Enable buffering. By default buffering for stderr and stdout is\n\
        disabled.\n\
\n\
    -v, --verbose\n\
        Be more verbose.\n\
\n\
    -h,--help\n\
        Show this help\n\
\n\
\n\
Compiled %s, %s\n", PRECISION, __DATE__, __TIME__);
}

int main(int argc, char *argv[])
{
    double precision = PRECISION;
    double T,Q;
    double lfac = 5;
    int i, n = -1, m = -1;
    int lmax = 0;
    int buffering_flag = 0, verbose_flag = 0;

    printf("# %s", argv[0]);
    for(i = 1; i < argc; i++)
        printf(", %s", argv[i]);
    printf("\n");

    while (1)
    {
        int c;
        struct option long_options[] =
        {
          { "verbose",   no_argument,       &verbose_flag,   1 },
          { "buffering", no_argument,       &buffering_flag, 1 },
          { "help",      no_argument,       0, 'h' },
          { "lscale",    required_argument, 0, 'l' },
          { "precision", required_argument, 0, 'p' },
          { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
      
        c = getopt_long (argc, argv, "Q:T:n:m:s:a:l:L:p:vqh", long_options, &option_index);
      
        /* Detect the end of the options. */
        if (c == -1)
          break;
      
        switch (c)
        {
          case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
              break;
          case 'Q':
              Q = atof(optarg);
              break;
          case 'T':
              T = atof(optarg);
              break;
          case 'L':
              lmax = atoi(optarg);
          case 'v':
              verbose_flag = 1;
              break;
          case 'l':
              lfac = atof(optarg);
              break;
          case 'p':
              precision = atof(optarg);
              break;
          case 'n':
              n = atoi(optarg);
              break;
          case 'm':
              m = atoi(optarg);
              break;
          case 'h':
              usage(stdout);
              exit(0);
      
          case '?':
            /* getopt_long already printed an error message. */
            break;
      
          default:
            abort();
        }
    }

    // disable buffering
    if(!buffering_flag)
    {
        fflush(stdin);
        fflush(stderr);
        setvbuf(stdout, NULL, _IONBF, 0);
        setvbuf(stderr, NULL, _IONBF, 0);
    }

    if(lfac <= 0)
    {
        fprintf(stderr, "--lfac must be positive\n\n");
        usage(stderr);
        exit(1);
    }
    if(precision <= 0)
    {
        fprintf(stderr, "--precision must be positive\n\n");
        usage(stderr);
        exit(1);
    }
    if(Q <= 0 || Q >= 1)
    {
        fprintf(stderr, "-Q must be in 0 < Q < 1\n\n");
        usage(stderr);
        exit(1);
    }
    if(T <= 0)
    {
        fprintf(stderr, "positive value for -T required\n\n");
        usage(stderr);
        exit(1);
    }
    if(n < 0 || m < 0)
    {
        fprintf(stderr, "n,m >= 0\n\n");
        usage(stderr);
        exit(1);
    }

    printf("# precision=%g\n", precision);
    printf("# lfac=%g\n", lfac);
    printf("# Q = %g\n", Q);
    printf("# n = %d\n", n);
    printf("# m = %d\n", m);
    printf("#\n");
    printf("# R/(L+R), T, F, lmax, nmax, time\n");

    {
        casimir_t casimir;
        casimir_mie_cache_t cache;
        double value, start_time = now();
        double nTRbyScriptL = n*T*Q;

        casimir_init(&casimir, Q, T);
        casimir_set_precision(&casimir, precision);
        casimir_set_verbose(&casimir, verbose_flag);

        if(lmax > 0)
            casimir_set_lmax(&casimir, lmax);
        else
            casimir_set_lmax(&casimir, MAX((int)ceil(Q/(1-Q)*lfac), 5));

        casimir_mie_cache_init(&cache, nTRbyScriptL);
        casimir_mie_cache_alloc(&casimir, &cache, casimir.lmax);
        value = casimir_logdetD(&casimir, n, m, &cache);
        casimir_mie_cache_free(&cache);
        casimir_free(&casimir);

        printf("# Q,T,n,m,value,lmax,time\n");
        printf("%.15g, %.15g, %d, %d, %.15g, %d, %g\n", Q, T, n, m, value, casimir.lmax, now()-start_time);
    }

    return 0;
}
