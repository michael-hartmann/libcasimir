#include <ctype.h>
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

#define PRECISION 1e-10

/* This function counts how many times the character c in string str appears */
int count(const char *str, char c)
{
    int i = 0;
    while(*str++ != '\0')
        if(*str == c)
            i++;

    return i;
}

/* This function returns a pointer to the n-th occurrence of the character c in
 * the string s. If the character c occures less than n times, NULL is
 * returned. 
 */
const char *indexn(const char *str, char c, int n)
{
    int i = 0;
    while(*str++ != '\0')
        if(*str == c && ++i == n)
            return str;

    return NULL;
}

/* This function returns the seconds since 1st Jan 1970 in µs precision */
double now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec + tv.tv_usec*1e-6;
}

/* print usage */
void usage(FILE *stream)
{
    fprintf(stream, "Usage: casimir [OPTIONS]\n\
This program will calculate the free Casimir energy for the plane-sphere \n\
geometry. \n\
\n\
Mandatory options:\n\
    -Q Radius of sphere divided by distance between center of sphere and plate; 0 < R/(R+L) < 1\n\
    -T Temperature\n\
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
    -q, --quiet\n\
        The progress is printed to stderr unless this flag is set.\n\
\n\
    -h,--help\n\
        Show this help\n\
\n\
\n\
Compiled %s, %s\n", PRECISION, __DATE__, __TIME__);
}

double iv(double list[4], int i)
{
    if(list[2] == 1)
        return list[0];

    if(list[3] == SCALE_LIN)
        return list[0]+(list[1]-list[0])*i/(list[2]-1);
    else
        return list[0]*pow(pow(list[1]/list[0], 1/(list[2]-1)), i);
}

void parse_range(char param, const char *_optarg, double list[])
{
    int elems = count(_optarg, ',');
    list[3] = SCALE_LIN;

    switch(elems)
    {
        case 0:
            list[0] = list[1] = atof(_optarg);
            list[2] = 1;
            break;
        case 3:
            if(strncasecmp(indexn(_optarg, ',', 3)+1, "log", 3) == 0)
                list[3] = SCALE_LOG;
        case 2:
            list[0] = atof(_optarg);
            list[1] = atof(indexn(_optarg, ',', 1)+1);
            list[2] = atoi(indexn(_optarg, ',', 2)+1);

            if(list[0] < 0 || list[1] < 0 || list[2] < 0)
            {
                fprintf(stderr, "error parsing parameter -%c\n\n", param);
                usage(stderr);
                exit(1);
            }

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
    double precision = PRECISION;
    double lT[4] = { 0,0,0,SCALE_LIN };
    double lQ[4] = { 0,0,0,SCALE_LIN };
    double lfac = 5;
    int i, iT, iQ;
    int lmax = 0;
    int buffering_flag = 0, quiet_flag = 0, verbose_flag = 0;

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
          { "quiet",     no_argument,       &quiet_flag,     1 },
          { "buffering", no_argument,       &buffering_flag, 1 },
          { "help",      no_argument,       0, 'h' },
          { "lscale",    required_argument, 0, 'l' },
          { "precision", required_argument, 0, 'p' },
          { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
      
        c = getopt_long (argc, argv, "Q:T:s:a:l:L:p:vqh", long_options, &option_index);
      
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
              parse_range('Q', optarg, lQ);
              break;
          case 'T':
              parse_range('T', optarg, lT);
              break;
          case 'L':
              lmax = atoi(optarg);
          case 'q':
              quiet_flag = 1;
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
    if(lQ[0] <= 0 || lQ[0] >= 1 || lQ[1] <= 0 || lQ[1] >= 1)
    {
        fprintf(stderr, "-Q must be in 0 < Q < 1\n\n");
        usage(stderr);
        exit(1);
    }
    if(lT[0] <= 0)
    {
        fprintf(stderr, "positive value for -T required\n\n");
        usage(stderr);
        exit(1);
    }

    printf("# precision=%g\n", precision);
    printf("# lfac=%g\n", lfac);
    if(lQ[2] == 1)
        printf("# L=%g\n", lQ[0]);
    else
        printf("# L=%g...%g (%d)\n", lQ[0],lQ[1],(int)lQ[2]);
    if(lT[2] == 1)
        printf("# T=%g\n", lT[0]);
    else
        printf("# T=%g...%g (%d)\n", lT[0],lT[1],(int)lT[2]);

    printf("#\n");
    printf("# R/(L+R), T, F, lmax, nmax, time\n");

    i = 0;
    for(iQ = 0; iQ < lQ[2]; iQ++)
        for(iT = 0; iT < lT[2]; iT++)
        {
            casimir_t casimir;
            double start_time = now();
            int nmax;
            double F,Q,T;

            Q = iv(lQ, iQ);
            T = iv(lT, iT);

            casimir_init(&casimir, Q, T);
            casimir_set_precision(&casimir, precision);
            casimir_set_verbose(&casimir, verbose_flag);

            if(lmax > 0)
                casimir_set_lmax(&casimir, lmax);
            else
                casimir_set_lmax(&casimir, MAX((int)ceil(Q/(1-Q)*lfac), 5));

            F = casimir_F(&casimir, &nmax);
            casimir_free(&casimir);

            printf("%.15g, %.15g, %.15g, %d, %d, %g\n", Q, T, F, casimir.lmax, nmax, now()-start_time);

            if(!quiet_flag)
                fprintf(stderr, "# %6.2f%%, R/(R+L)=%g, T=%g\n", ++i*100/(lQ[2]*lT[2]), Q, T);
        }

    return 0;
}
