#ifndef __CASIMIR_H
#define __CASIMIR_H

    #define SCALE_LIN 0
    #define SCALE_LOG 1

    /* prototypes */
    int count(const char *str, char c);
    const char *indexn(const char *str, char c, int n);
    double now(void);
    void usage(FILE *stream);
    double iv(double list[4], int i);
    void parse_range(char param, const char *_optarg, double list[]);

#endif
