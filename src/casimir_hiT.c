#define _GNU_SOURCE

#include <getopt.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "casimir_hiT.h"
#include "libcasimir.h"

#define IDLE 1000 // in µs
#define LSCALE 7.0
#define PRECISION 5e-9
#define VERBOSE 1

static void *xmalloc(size_t size)
{
    void *p = malloc(size);
    if(p == NULL)
    {
        fprintf(stderr, "Out of memory: malloc failed.\n");
        exit(2);
    }
    return p;
}

/* This function returns the seconds since 1st Jan 1970 in µs precision */
double now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec + tv.tv_usec*1e-6;
}

double sumF(double *values, int lmax)
{
    int i;
    double F = 0;

    for(i = lmax-1; i >= 0; i--)
        F += values[i];

    return F;
}

void usage(FILE *stream, const char *name)
{
    fprintf(stream, "Usage: %s -x L/R [-l lscale -L lmax -p precision]\n\n", name);
    fprintf(stream, "\t-x L/R:    ratio of L and R, L/R > 0\n");
    fprintf(stream, "\t-L lmax:   use lmax\n"); 
    fprintf(stream, "\t-l lscale: use lmax = lscale*R/L (ignored if -L is used), default: %g\n", LSCALE); 
    fprintf(stream, "\t-p prec:   use precision, default: %g\n", PRECISION); 
    fprintf(stream, "\t-c cores:  how many cores to use, default: 1\n");
}

pthread_t *start_thread(double LbyR, int m, int lmax, double precision)
{
    pthread_t *thread = xmalloc(sizeof(pthread_t));
    param_t *p        = xmalloc(sizeof(param_t));

    p->LbyR      = LbyR;
    p->precision = precision;
    p->lmax      = lmax;
    p->m         = m;
    p->value     = -1;
    p->time      = -1;

    pthread_create(thread, NULL, &logdetD0, (void *)p);

    return thread;
}


void *logdetD0(void *p)
{       
    casimir_t casimir;
    double start = now();
    double value;
    param_t *params = p;
    double LbyR      = params->LbyR;
    double precision = params->precision;
    int m            = params->m;
    int lmax         = params->lmax;

    casimir_init(&casimir, 1/(1+LbyR), 0.1);
    casimir_set_precision(&casimir, precision);
    casimir_set_verbose(&casimir, VERBOSE);
    casimir_set_lmax(&casimir, lmax);
    
    value = casimir_logdetD(&casimir,0,m,NULL)/2;
    casimir_free(&casimir);
    
    if(m == 0)
        value /= 2;
    
    params->value = value;
    params->time  = now()-start;

    return params;
}

int main(int argc, char *argv[])
{
    double precision = PRECISION;
    double lscale = LSCALE;
    int lmax = -1;
    int cores = 1;
    int i, m;
    double start_time = now();
    double LbyR = -1;
    double *values;
    pthread_t **threads;

    while (1)
    {
        int c;
        struct option long_options[] =
        {
          { "help",      no_argument,       0, 'h' },
          { "LbyR",      required_argument, 0, 'x' },
          { "lmax",      required_argument, 0, 'L' },
          { "lscale",    required_argument, 0, 'l' },
          { "precision", required_argument, 0, 'p' },
          { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
      
        c = getopt_long (argc, argv, "x:L:l:p:c:h", long_options, &option_index);
      
        /* Detect the end of the options. */
        if (c == -1)
          break;
      
        switch (c)
        {
          case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
              break;
          case 'x':
              LbyR = atof(optarg);
              break;
          case 'L':
              lmax = atoi(optarg);
              break;
          case 'l':
              lscale = atof(optarg);
              break;
          case 'c':
              cores = atoi(optarg);
              break;
          case 'p':
              precision = atof(optarg);
              break;
          case 'h':
              usage(stdout, argv[0]);
              exit(0);
      
          case '?':
            /* getopt_long already printed an error message. */
            break;
      
          default:
            abort();
        }
    }

    if(LbyR <= 0)
    {
        fprintf(stderr, "argument of -x must be nonnegative\n\n");
        usage(stderr, argv[0]);
        exit(1);
    }
    if(precision <= 0)
    {
        fprintf(stderr, "argument of -p must be nonnegative\n\n");
        usage(stderr, argv[0]);
        exit(1);
    }
    if(cores < 0)
    {
        fprintf(stderr, "argument of -t must be at least 1\n\n");
        usage(stderr, argv[0]);
        exit(1);
    }
    if(lmax <= 0)
        lmax = lscale/LbyR;

    // disable buffering
    {
        fflush(stdin);
        fflush(stderr);
        setvbuf(stdout, NULL, _IONBF, 0);
        setvbuf(stderr, NULL, _IONBF, 0);
    }

    fprintf(stderr, "# L/R=%g, precision=%g, lmax=%d, cores=%d\n", LbyR, precision, lmax, cores);
    
    values = (double *)xmalloc(lmax*sizeof(double));

    for(i = 0; i < lmax; i++)
        values[i] = 0;

    threads = (pthread_t **)xmalloc(cores*sizeof(pthread_t));
    for(i = 0; i < cores; i++)
        threads[i] = NULL;

    m = 0;
    while(m < lmax)
    {
        int bye = 0;
        param_t *r;

        // try to start threads
        for(i = 0; i < cores; i++)
        {
            if(threads[i] == NULL)
                threads[i] = start_thread(LbyR, m++, lmax, precision);
            else if(pthread_tryjoin_np(*threads[i], (void *)&r) == 0)
            {
                values[r->m] = r->value;
                fprintf(stderr, "# m=%d, value=%.15g, time=%g\n", r->m, r->value, r->time);
                if(r->value/sumF(values, lmax) < precision)
                    bye = 1;
                free(r);
                threads[i] = NULL;
            }
        }

        if(bye)
            break;

        usleep(IDLE);
    }

    fprintf(stderr, "# waiting for last threads to finish\n");
    for(i = 0; i < cores; i++)
    {
        param_t *r;
        if(threads[i] != NULL)
        {
            pthread_join(*threads[i], (void *)&r);
            values[r->m] = r->value;
            fprintf(stderr, "# m=%d, value=%.15g, time=%g\n", r->m, r->value, r->time);
            free(r);
            threads[i] = NULL;
        }
    }
    
    free(threads);

    printf("# L/R, value, time\n");
    printf("%.15g, %.15g, %.15g\n", LbyR, sumF(values, lmax), now()-start_time);

    free(values);

    return 0;
}
