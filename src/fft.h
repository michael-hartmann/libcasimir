#ifndef FFT_H
#define FFT_H

typedef struct {
    int n;
    int log2n;
    double *W_re;
    double *W_im;
} fft_t;

void fft_init(fft_t *self, int n);
void fft_free(fft_t *self);

void fft(fft_t *self, int dir, int n, double *A_re, double *A_im);

int  fft_bitrev(int inp, int numbits); 
void fft_permute_bitrev(fft_t *self, double *A_re, double *A_im);
void fft_compute_W(int n, double *W_re, double *W_im); ;

#endif
