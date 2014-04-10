#ifndef UNIFFT_H
#define UNIFFT_H

#include "edouble.h"

typedef struct {
    int n, log2n;
    edouble *W_re, *W_im;
} fft_t;

void fft_permute_bitrev(fft_t *self, edouble *A_re, edouble *A_im);
int  fft_fft(fft_t *self, int dir, edouble *A_re, edouble *A_im);
int  fft_init(fft_t *self, int n);
void fft_free(fft_t *self);
int fft_polymult(edouble p1[], edouble p2[], size_t len);

#endif
