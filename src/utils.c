#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include "utils.h"

void (*error_handler)(const char *) = default_error_handler;

void set_defaut_error_handler(void (*f)(const char *))
{
    error_handler = f;
}

void default_error_handler(const char *reason)
{
    fprintf(stderr, "%s\n", reason);
    exit(1);
}

void *xmalloc(size_t len)
{
    void *p = malloc(len);
    if(p == NULL)
        error_handler("malloc failed.");

    return p;
}

void *xrealloc(void *p, size_t len)
{
    void *p2 = realloc(p, len);
    if(p2 == NULL)
    {
        free(p);
        error_handler("realloc failed.");
    }

    return p2;
}

void xfree(void *p)
{
    if(p == NULL)
        error_handler("free on NULL.");

    free(p);
}

/* This function returns the seconds since 1st Jan 1970 in µs precision */
double now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec + tv.tv_usec*1e-6;
}

