#include <stdio.h>
#include <ctype.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include "libcasimir.h"

#define SCALE_LIN 0
#define SCALE_LOG 1

#define kB      1.38064880e-23
#define HBAR    1.05457173e-34
#define C       299792458

int count(const char *str, char c)
{
    int i = 0;
    while(*str++ != '\0')
        if(*str == c)
            i++;

    return i;
}

const char *indexn(const char *str, char c, int n)
{
    int i = 0;
    while(*str++ != '\0')
        if(*str == c)
            if(++i == n)
                return str;

    return NULL;
}

double now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec + tv.tv_usec*1e-6;
}

void usage(FILE *stream)
{
    fprintf(stream, "Usage: casimir [OPTIONS]\n\
This program will calculate the free Casimir energy for the plane-sphere \n\
geometry. \n\
\n\
Mandatory options:\n\
    -R  Radius of sphere\n\
    -L  Distance between sphere and plate\n\
    -T  Temperature\n\
\n\
Further options:\n\
    -l, --lscale\n\
        Specify parameter lscale. The vector space has to be truncated at some\n\
        value lmax. This program will use lmax=(R/L*lscale) (default: 3)\n\
\n\
    --buffering\n\
        Enable buffering. By default buffering for stderr and stdout is\n\
        disabled.\n\
\n\
    --quiet\n\
        The progress and the ETA is printed to stderr unless this flag is set.\n\
\n\
    -h,--help\n\
        Show this help\n");
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

void parse_range(char param, const char *optarg, double list[])
{
    int elems = count(optarg, ',');
    list[3] = SCALE_LIN;

    switch(elems)
    {
        case 0:
            list[0] = list[1] = atof(optarg);
            list[2] = 1;
            break;
        case 3:
            if(strncasecmp(indexn(optarg, ',', 3)+1, "log", 3) == 0)
                list[3] = SCALE_LOG;
        case 2:
            list[0] = atof(optarg);
            list[1] = atof(indexn(optarg, ',', 1)+1);
            list[2] = atoi(indexn(optarg, ',', 2)+1);

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
            fprintf(stderr, "Can't parse range %s.\n\n", optarg);
            usage(stderr);
            exit(1);
    }
}

int main(int argc, char *argv[])
{
    double lT[4] = { 0,0,0,SCALE_LIN };
    double lR[4] = { 0,0,0,SCALE_LIN };
    double lL[4] = { 0,0,0,SCALE_LIN };
    double lfac = 3;
    int i, iR, iL, iT;
    double R,L,T;
    int lmax;
    int buffering_flag = 0, quiet_flag = 0;

    printf("# %s", argv[0]);
    for(i = 1; i < argc; i++)
        printf(", %s", argv[i]);
    printf("\n");

    while (1)
    {
        int c;
        struct option long_options[] =
        {
          { "quiet",     no_argument,       &quiet_flag,     1 },
          { "buffering", no_argument,       &buffering_flag, 1 },
          { "help",      no_argument,       0, 'h' },
          { "lscale",    required_argument, 0, 'l' },
          { 0, 0, 0, 0 }
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
      
        c = getopt_long (argc, argv, "R:L:T:a:l:h", long_options, &option_index);
      
        /* Detect the end of the options. */
        if (c == -1)
          break;
      
        switch (c)
        {
          case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
              break;
            printf ("option %s", long_options[option_index].name);
            if (optarg)
              printf (" with arg %s", optarg);
            printf ("\n");
            break;

          case 'R':
              parse_range('R', optarg, lR);
              break;
          case 'L':
              parse_range('L', optarg, lL);
              break;
          case 'T':
              parse_range('T', optarg, lT);
              break;
          case 'l':
              lfac = atof(optarg);
              break;
          case 'h':
              usage(stdout);
              exit(0);
      
          case '?':
            /* getopt_long already printed an error message. */
            break;
      
          default:
            abort ();
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
    if(lR[0] <= 0)
    {
        fprintf(stderr, "positive value for -R required\n\n");
        usage(stderr);
        exit(1);
    }
    if(lL[0] <= 0)
    {
        fprintf(stderr, "positive value for -L required\n\n");
        usage(stderr);
        exit(1);
    }
    if(lT[0] <= 0)
    {
        fprintf(stderr, "positive value for -T required\n\n");
        usage(stderr);
        exit(1);
    }

    printf("# lfac=%g, ", lfac);
    if(lR[2] == 1)
        printf("R=%g, ", lR[0]);
    else
        printf("R=%g...%g (%d), ", lR[0],lR[1],(int)lR[2]);
    if(lL[2] == 1)
        printf("L=%g, ", lL[0]);
    else
        printf("L=%g...%g (%d), ", lL[0],lL[1],(int)lL[2]);
    if(lT[2] == 1)
        printf("T=%g\n", lT[0]);
    else
        printf("T=%g...%g (%d)\n", lT[0],lT[1],(int)lT[2]);

    printf("# R, L, T, F, L/R, lmax, nmax, time\n");

    i = 0;
    for(iL = 0; iL < lL[2]; iL++)
        for(iR = 0; iR < lR[2]; iR++)
            for(iT = 0; iT < lT[2]; iT++)
            {
                casimir_t casimir;
                double start_time = now();
                double F, F_scaled, T_scaled;
                int nmax;

                R = iv(lR, iR);
                L = iv(lL, iL);
                T = iv(lT, iT);

                T_scaled = T*(L+R) * (2*M_PI*kB)/(HBAR*C);

                lmax = MAX((int)ceil(R/L*lfac), 3);

                casimir_init_perfect(&casimir, R/(R+L), T_scaled);
                casimir_set_lmax(&casimir, lmax);
                F_scaled = casimir_F(&casimir, &nmax);
                F = F_scaled*HBAR*C/(L+R);
                casimir_free(&casimir);

                printf("%.15g, %.15g, %.15g, %.15g, %.15g, %d, %d, %g\n", R, L, T, F, L/R, lmax, nmax, now()-start_time);

                if(!quiet_flag)
                    fprintf(stderr, "# %6.2f%%, R=%g, L=%g, T=%g\n", ++i*100/(lL[2]*lR[2]*lT[2]), R, L, T);
            }

    return 0;
}
