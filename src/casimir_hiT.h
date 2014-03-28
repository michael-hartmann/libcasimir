#ifndef CASIMIR_HIT
#define CASIMIR_HIT

typedef struct {
    int m;
    double value;
    double time;
} return_t;

typedef struct {
    double LbyR;
    int m;
    double precision;
    int lmax;
} param_t;

double now(void);
double sumF(double *values, int lmax);
void usage(FILE *stream, const char *name);
void *logdetD0(void *p);
pthread_t *start_thread(double LbyR, int m, int lmax, double precision);

#endif
