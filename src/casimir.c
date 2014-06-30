#include <ctype.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "casimir.h"
#include "edouble.h"
#include "libcasimir.h"
#include "sfunc.h"
#include "utils.h"

/* default values for precision and lfac */
#define DEFAULT_PRECISION 1e-10
#define DEFAULT_LFAC 5


/* print usage */
void usage(FILE *stream)
{
    fprintf(stream, "Usage: casimir [OPTIONS]\n\
This program will calculate the free Casimir energy F(T,L/R) for the\n\
plane-sphere geometry for given L/R and temperature T. The output is in scaled\n\
unities.\n\
\n\
Mandatory options:\n\
    -x, --LbyR FRACTION\n\
        Separation L between sphere and plane divided by radius of of sphere,\n\
        where L/R > 0.\n\
        If you want to calculate several points, you may pass start and stop\n\
        value and the amount of points to be calculated.\n\
        Examples:\n\
            $ ./casimir -T 1 -x 0.5,0.9,5\n\
            This will calculate five free energies for Q=0.5,...,0,9 in linear\n\
            scale.\n\
            $ ./casimir -T 1 -x 0.5,0.9,5,log\n\
            This will calculate five free energies for Q=0.5,...,0,9, but using\n\
            a logarithmic scale.\n\
\n\
    -T TEMPERATURE\n\
        Temperature in scaled units. You may use the same syntax as for -Q.\n\
\n\
Further options:\n\
    -l, --lscale\n\
        Specify parameter lscale. The vector space has to be truncated for\n\
        some value lmax. This program will use lmax=(R/L*lscale) (default: %d)\n\
\n\
    -L LMAX\n\
        Set lmax to the value LMAx. When -L is specified, -l will be ignored\n\
\n\
    -c, --cores CORES\n\
        Use CORES at the same time (default: 1)\n\
\n\
    -p, --precision\n\
        Set precision to given value. (default: %g)\n\
\n\
    --buffering\n\
        Enable buffering. By default buffering for stderr and stdout is\n\
        disabled.\n\
\n\
    -v, --verbose\n\
        Be more verbose. This will output additional information.\n\
\n\
    -X, --extrapolate\n\
        Achieve better results by extrapolating the contributions F_n.\n\
        This feature is experimental! Use it on your own risk! This feature might\n\
        cause your computer to explode or even worse: It may cause wrong results!\n\
        Don't blame me, I have warned you!\n\
\n\
    -q, --quiet\n\
        The progress is printed to stderr unless this flag is set.\n\
\n\
    -h,--help\n\
        Show this help\n\
\n\
\n\
Compiled %s, %s\n\
%s\n", DEFAULT_LFAC, DEFAULT_PRECISION, __DATE__, __TIME__, casimir_compile_info());
}

/* parse a range given for LbyR or T from the command line.
 * Examples:
 * 1) "value"            => list = { value, value, 1, SCALE_LIN }
 * 2) "start,stop,N"     => list = { start, stop, N, SCALE_LIN }
 * 3) "start,stop,N,log" => list = { start, stop, N, SCALE_LOG }
 */
void parse_range(const char param, const char *_optarg, double list[])
{
    int elems = cinstr(_optarg, ','); /* commas in _optarg */
    list[3] = SCALE_LIN;

    switch(elems)
    {
        case 0:
            /* no comma => example 1) */
            list[0] = list[1] = atof(_optarg);
            list[2] = 1;
            break;
        case 3:
            /* 3 commas => example 3) */
            if(strncasecmp(indexn(_optarg, ',', 3)+1, "log", 3) == 0)
                list[3] = SCALE_LOG;
            /* here no break! */
        case 2:
            /* 2 commas => example 2) */
            list[0] = atof(_optarg);
            list[1] = atof(indexn(_optarg, ',', 1)+1);
            list[2] = atoi(indexn(_optarg, ',', 2)+1);

            if(list[0] < 0 || list[1] < 0 || list[2] < 0)
            {
                fprintf(stderr, "error parsing parameter -%c\n\n", param);
                usage(stderr);
                exit(1);
            }

            /* ensure that start < stop */
            if(list[0] > list[1])
            {
                double temp = list[0];
                list[0] = list[1];
                list[1] = temp;
            }
            break;

        default:
            fprintf(stderr, "Can't parse range %s.\n\n", _optarg);
            usage(stderr);
            exit(1);
    }
}


int main(int argc, char *argv[])
{
    double precision = DEFAULT_PRECISION;
    double lT[4] = { 0,0,0,SCALE_LIN }; /* start, stop, N, lin/log */
    double lLbyR[4] = { 0,0,0,SCALE_LIN }; /* start, stop, N, lin/log */
    double lfac = DEFAULT_LFAC;
    int i, iT, iLbyR;
    int cores = 1;
    int lmax = 0;
    int buffering_flag = 0, quiet_flag = 0, verbose_flag = 0, extrapolate_flag = 0;

    /* parse command line options */
    while (1)
    {
        int c;
        struct option long_options[] =
        {
          { "verbose",     no_argument,       &verbose_flag,     1 },
          { "quiet",       no_argument,       &quiet_flag,       1 },
          { "buffering",   no_argument,       &buffering_flag,   1 },
          { "extrapolate", no_argument,       &extrapolate_flag, 1 },
          { "help",        no_argument,       0, 'h' },
          { "LbyR",        required_argument, 0, 'x' },
          { "lscale",      required_argument, 0, 'l' },
          { "cores",       required_argument, 0, 'c' },
          { "precision",   required_argument, 0, 'p' },
          { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
      
        c = getopt_long (argc, argv, "x:T:c:s:a:l:L:p:Xvqh", long_options, &option_index);
      
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
                parse_range('x', optarg, lLbyR);
                break;
            case 'T':
                parse_range('T', optarg, lT);
                break;
            case 'L':
                lmax = atoi(optarg);
                break;
            case 'X':
                extrapolate_flag = 1;
                break;
            case 'q':
                quiet_flag = 1;
                break;
            case 'c':
                cores = atoi(optarg);
                break;
            case 'v':
                verbose_flag = 1;
                break;
            case 'l':
                lfac = atof(optarg);
                break;
            case 'p':
                precision = atof(optarg);
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

    /* check command line arguments */
    {
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
        if(lLbyR[0] <= 0 || lLbyR[1] <= 0)
        {
            fprintf(stderr, "-x must be positive; LbyR > 0\n\n");
            usage(stderr);
            exit(1);
        }
        if(lT[0] <= 0)
        {
            fprintf(stderr, "positive value for -T required\n\n");
            usage(stderr);
            exit(1);
        }
    }

    /* print information to stdout */
    {
        /* print command line */
        printf("# %s", argv[0]);
        for(i = 1; i < argc; i++)
            printf(", %s", argv[i]);
        printf("\n");

        printf("# precision=%g\n", precision);
        if(lmax > 0)
            printf("# lmax=%d\n", lmax);
        else
            printf("# lfac=%g\n", lfac);
        printf("# cores=%d\n", cores);
        if(lLbyR[2] == 1)
            printf("# LbyR=%g\n", lLbyR[0]);
        else
            printf("# LbyR=%g...%g (%d)\n", lLbyR[0],lLbyR[1],(int)lLbyR[2]);
        if(lT[2] == 1)
            printf("# T=%g\n", lT[0]);
        else
            printf("# T=%g...%g (%d)\n", lT[0],lT[1],(int)lT[2]);
        printf("# extrapolate=%s\n", extrapolate_flag ? "yes" : "no");

        printf("#\n");
        printf("# LbyR, T, F, lmax, nmax, time\n");
    }

    i = 0;
    for(iLbyR = 0; iLbyR < lLbyR[2]; iLbyR++)
        for(iT = 0; iT < lT[2]; iT++)
        {
            casimir_t casimir;
            double start_time = now();
            int nmax;
            double F,LbyR,T,Q;

            if(lLbyR[3] == SCALE_LIN)
                LbyR = linspace(lLbyR[0], lLbyR[1], lLbyR[2], iLbyR);
            else
                LbyR = logspace(lLbyR[0], lLbyR[1], lLbyR[2], iLbyR);

            if(lT[3] == SCALE_LIN)
                T = linspace(lT[0], lT[1], lT[2], iT);
            else
                T = logspace(lT[0], lT[1], lT[2], iT);

            Q = 1/(1+LbyR);
            casimir_init(&casimir, Q, T);
            casimir_set_cores(&casimir, cores);
            casimir_set_precision(&casimir, precision);
            casimir_set_verbose(&casimir, verbose_flag);
            casimir_set_extrapolate(&casimir, extrapolate_flag);

            if(lmax > 0)
                casimir_set_lmax(&casimir, lmax);
            else
                casimir_set_lmax(&casimir, MAX((int)ceil(lfac/LbyR), DEFAULT_LFAC));

            F = casimir_F(&casimir, &nmax);
            casimir_free(&casimir);

            printf("%.15g, %.15g, %.15g, %d, %d, %g\n", LbyR, T, F, casimir.lmax, nmax, now()-start_time);

            if(!quiet_flag)
                fprintf(stderr, "# %6.2f%%, L/R=%g, T=%g\n", ++i*100/(lLbyR[2]*lT[2]), LbyR, T);
        }

    return 0;
}
