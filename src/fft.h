#ifndef FFT_H
#define FFT_H

#include "integration.h"

typedef struct {
    int n, log2n;
    log_t *W_re, *W_im;
} fft_t;

int  fft_init(fft_t *self, int n);
void fft_free(fft_t *self);

void fft(fft_t *self, int dir, log_t *A_re, log_t *A_im);
void fft_permute_bitrev(fft_t *self, log_t *A_re, log_t *A_im);

#endif
