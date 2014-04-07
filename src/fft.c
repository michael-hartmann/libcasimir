#include <stdlib.h> 
#include <math.h>   
#include <stdio.h>

#include "fft.h"

/* returns log n (to the base 2), if n is positive and power of 2 */ 
static int fft_log_2(int n) 
{
    int res; 
    for(res=0; n >= 2; res++) 
        n = n >> 1; 

    return res; 
}
 

/* treats inp as a numbits number and bitreverses it. 
 * inp < 2^(numbits) for meaningful bit-reversal
 */ 
static int fft_bitrev(int inp, int numbits)
{
    int i, rev=0;

    for(i = 0; i < numbits; i++)
    {
        rev = (rev << 1) | (inp & 1);
        inp >>= 1;
    }

    return rev;
}

/* permutes the array using a bit-reversal permutation */ 
void fft_permute_bitrev(fft_t *self, double *A_re, double *A_im) 
{ 
    int i, bri;
    double t_re, t_im;
    int n     = self->n;
    int log2n = self->log2n;
    
    for(i = 0; i < n; i++)
    {
        bri = fft_bitrev(i, log2n);
    
        /* skip already swapped elements */
        if (bri <= i)
            continue;
    
        t_re      = A_re[i];
        t_im      = A_im[i];
        A_re[i]   = A_re[bri];
        A_im[i]   = A_im[bri];
        A_re[bri] = t_re;
        A_im[bri] = t_im;
    }
} 


/* W will contain roots of unity so that W[bitrev(i,log2n-1)] = e^(2*pi*i/n)
 * n should be a power of 2
 * Note: W is bit-reversal permuted because fft(..) goes faster if this is done.
 *       see that function for more details on why we treat 'i' as a (log2n-1) bit number.
 */
void fft_init(fft_t *self, int n)
{
    int i, br;
    int log2n = fft_log_2(n);

    self->n     = n;
    self->log2n = log2n;

    self->W_re = (double *)malloc(n/2*sizeof(double));
    self->W_im = (double *)malloc(n/2*sizeof(double));
    
    for (i = 0; i < (n/2); i++)
    {
        br = fft_bitrev(i, log2n-1); 
        self->W_re[br] = +cos(i*2*M_PI/n);  
        self->W_im[br] = -sin(i*2*M_PI/n);  
    }
}

void fft_free(fft_t *self)
{
    if(self->W_re != NULL)
        free(self->W_re);
    if(self->W_im != NULL)
        free(self->W_im);

    self->W_re = self->W_im = NULL;
}




/* fft on a set of n points given by A_re and A_im. Bit-reversal permuted roots-of-unity lookup table
 * is given by W_re and W_im. More specifically,  W is the array of first n/2 nth roots of unity stored
 * in a permuted bitreversal order.
 *
 * FFT - Decimation In Time FFT with input array in correct order and output array in bit-reversed order.
 *
 * REQ: n should be a power of 2 to work. 
 *
 * Note: - See www.cs.berkeley.edu/~randit for her thesis on VIRAM FFTs and other details about VHALF section of the algo
 *         (thesis link - http://www.cs.berkeley.edu/~randit/papers/csd-00-1106.pdf)
 *       - See the foll. CS267 website for details of the Decimation In Time FFT implemented here.
 *         (www.cs.berkeley.edu/~demmel/cs267/lecture24/lecture24.html)
 *       - Also, look "Cormen Leicester Rivest [CLR] - Introduction to Algorithms" book for another variant of Iterative-FFT
 */

void fft(fft_t *self, int dir, int n, double *A_re, double *A_im) 
{
    double w_re, w_im, u_re, u_im, t_re, t_im;
    int m, g, b;
    int mt, k;
    double *W_re = self->W_re;
    double *W_im = self->W_im;
    
    if(dir > 0)
        dir = 1;
    else
        dir = -1;
    
    /* for each stage */  
    for(m = n; m >= 2; m = m>>1)
    {
      /* m = n/2^s; mt = m/2; */
      mt = m >> 1;
    
      /* for each group of butterfly */ 
      for (g = 0, k = 0; g < n; g += m, k++) 
      {
        /* each butterfly group uses only one root of unity. actually, it is the bitrev of this group's number k.
         * BUT 'bitrev' it as a log2n-1 bit number because we are using a lookup array of nth root of unity and
         * using cancellation lemma to scale nth root to n/2, n/4,... th root.
         *
         * It turns out like the foll.
         *   w.re = W[bitrev(k, log2n-1)].re;
         *   w.im = W[bitrev(k, log2n-1)].im;
         * Still, we just use k, because the lookup array itself is bit-reversal permuted. 
         */
        w_re = W_re[k];
        w_im = W_im[k]*dir;
    
        /* for each butterfly */ 
        for (b=g; b<(g+mt); b++) 
        {
            /* t = w * A[b+mt] */
            t_re = w_re * A_re[b+mt] - w_im * A_im[b+mt];
            t_im = w_re * A_im[b+mt] + w_im * A_re[b+mt];
    
            /* u = A[b]; in[b] = u + t; in[b+mt] = u - t; */
            u_re = A_re[b];
            u_im = A_im[b];
            A_re[b] = u_re + t_re;
            A_im[b] = u_im + t_im;
            A_re[b+mt] = u_re - t_re;
            A_im[b+mt] = u_im - t_im;
        }
      }
    }
    
    if(dir == -1)
      for(m = 0; m < n; m++)
      {
          A_re[m] /= n;
          A_im[m] /= n;
      }
    
}

int main(int argc, char *argv[])
{
    int n;
    int i;
    fft_t self;
    
    n = 8;
    
    double A_re[]  = { 1, 2, 3, 4, 5, 6, 7, 8 };
    double A_im[]  = { 0, 0, 0, 0, 0, 0, 0, 0 };
    
    fft_init(&self, n);
    fft(&self, -1, n, A_re, A_im);
    fft_permute_bitrev(&self, A_re, A_im);        
    fft_free(&self);
    
    for(i = 0; i < n; i++)
        printf("%d: %+g%+gi\n", i, A_re[i], A_im[i]);
      
    return 0;
}
