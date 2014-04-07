#include <stdlib.h> 
#include <math.h>   
#include <stdio.h>

#include "fft.h"
#include "sfunc.h"
#include "integration.h"

/* returns log n (to the base 2), if n is positive and power of 2 */ 
static int fft_log_2(int n) 
{
    int res; 
    for(res = 0; n >= 2; res++) 
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
void fft_permute_bitrev(fft_t *self, log_t *A_re, log_t *A_im)
{
    int i;
    int n     = self->n;
    int log2n = self->log2n;
    
    for(i = 0; i < n; i++)
    {
        double t_re, t_im;
        int t_re_sign, t_im_sign; 
        int bri = fft_bitrev(i, log2n);
    
        /* skip already swapped elements */
        if(bri <= i)
            continue;
    
        t_re            = A_re[i].value;
        t_im            = A_im[i].value;
        A_re[i].value   = A_re[bri].value;
        A_im[i].value   = A_im[bri].value;
        A_re[bri].value = t_re;
        A_im[bri].value = t_im;

        t_re_sign = A_re[i].sign;
        t_im_sign = A_im[i].sign;
        A_re[i].sign = A_re[bri].sign;
        A_im[i].sign = A_im[bri].sign;
        A_re[bri].sign = t_re_sign;
        A_im[bri].sign = t_im_sign;
    }
} 


/* W will contain roots of unity so that W[bitrev(i,log2n-1)] = e^(2*pi*i/n)
 * n should be a power of 2
 * Note: W is bit-reversal permuted because fft(..) goes faster if this is done.
 *       see that function for more details on why we treat 'i' as a (log2n-1) bit number.
 */
int fft_init(fft_t *self, int n)
{
    int i, br;
    int log2n = fft_log_2(n);

    /* check if n is power of 2 */
    if(n != pow(2,log2n))
    {
        /* if it fails, at least set pointer to null, so we don't overwrite
         * other memory when the user doesn't check this return value */
        self->W_re = self->W_im = NULL;
        return 1;
    }

    self->n     = n;
    self->log2n = log2n;

    self->W_re = (log_t *)malloc(n/2*sizeof(log_t));
    self->W_im = (log_t *)malloc(n/2*sizeof(log_t));
    
    for (i = 0; i < (n/2); i++)
    {
        double re,im;
        br = fft_bitrev(i, log2n-1); 
        re = +cos(i*2*M_PI/n);
        im = -sin(i*2*M_PI/n);

        self->W_re[br].value = log(fabs(re));
        self->W_im[br].value = log(fabs(im));
        self->W_re[br].sign  = copysign(1, re);
        self->W_im[br].sign  = copysign(1, im);
    }

    return 0;
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
void fft(fft_t *self, int dir, log_t *A_re, log_t *A_im) 
{
    int n = self->n;
    int m, g, b, mt, k;
    
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
            double w_re = self->W_re[k].value;
            double w_im = self->W_im[k].value;
            int w_re_sign = self->W_re[k].sign;
            int w_im_sign = self->W_im[k].sign*dir;
            
            /* for each butterfly */ 
            for(b = g; b < (g+mt); b++) 
            {
                double u_re, u_im, t_re, t_im;
                int u_re_sign, u_im_sign, t_re_sign, t_im_sign;
                /* t = w * A[b+mt] */
                t_re = logadd_s(w_re+A_re[b+mt].value, w_re_sign*A_re[b+mt].sign, w_im+A_im[b+mt].value, -w_im_sign*A_im[b+mt].sign, &t_re_sign);
                t_im = logadd_s(w_re+A_im[b+mt].value, w_re_sign*A_im[b+mt].sign, w_im+A_re[b+mt].value, +w_im_sign*A_re[b+mt].sign, &t_im_sign);
            
                /* u = A[b]; in[b] = u + t; in[b+mt] = u - t; */
                u_re = A_re[b].value;
                u_im = A_im[b].value;
                u_re_sign = A_re[b].sign;
                u_im_sign = A_im[b].sign;

                A_re[b].value = logadd_s(u_re, u_re_sign, t_re, t_re_sign, &A_re[b].sign);
                A_im[b].value = logadd_s(u_im, u_im_sign, t_im, t_im_sign, &A_im[b].sign);

                A_re[b+mt].value = logadd_s(u_re, u_re_sign, t_re, -t_re_sign, &A_re[b+mt].sign);
                A_im[b+mt].value = logadd_s(u_im, u_re_sign, t_im, -t_im_sign, &A_im[b+mt].sign);
            }
        }
    }
    
    if(dir == -1)
    {
        double logn = log(n);
        for(m = 0; m < n; m++)
        {
            A_re[m].value -= logn;
            A_im[m].value -= logn;
        }
    }
}

int main(int argc, char *argv[])
{
    int n = 4;
    int i;
    fft_t self;
    log_t A_re[4];
    log_t A_im[4];
    
    A_re[0].value = log(1); A_re[0].sign = 1;
    A_re[1].value = log(2); A_re[1].sign = 1;
    A_re[2].value = log(3); A_re[2].sign = 1;
    A_re[3].value = log(4); A_re[3].sign = 1;

    A_im[0].value = log(0); A_im[0].sign = 1;
    A_im[1].value = log(0); A_im[1].sign = 1;
    A_im[2].value = log(0); A_im[2].sign = 1;
    A_im[3].value = log(0); A_im[3].sign = 1;
    
    fft_init(&self, n);
    fft(&self, -1, A_re, A_im);
    fft_permute_bitrev(&self, A_re, A_im);
    fft_free(&self);
    
    for(i = 0; i < n; i++)
        printf("%d: %+g%+gi\n", i, A_re[i].sign*exp(A_re[i].value), A_im[i].sign*exp(A_im[i].value));

    return 0;
}
