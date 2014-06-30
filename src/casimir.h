#ifndef __CASIMIR_H
#define __CASIMIR_H

    #define SCALE_LIN 0
    #define SCALE_LOG 1

    /* prototypes */
    void usage(FILE *stream);
    void parse_range(const char param, const char *_optarg, double list[]);

#endif
