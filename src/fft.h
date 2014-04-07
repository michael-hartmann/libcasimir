#ifndef FFT_H
#define FFT_H

void fft(int dir, int n, double *A_re, double *A_im, double *W_re, double *W_im);

int  fft_log_2(int n);  
int  fft_bitrev(int inp, int numbits); 
void fft_permute_bitrev(int n, double *A_re, double *A_im); 
void fft_compute_W(int n, double *W_re, double *W_im); ;

#endif
