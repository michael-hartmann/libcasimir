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

int fft_polymult3(log_t p1[], log_t p2[], log_t p3[], int len)
{
    int i;
    fft_t self;
    log_t *p1_im, *p2_im, *p3_im;

    if(fft_init(&self, len) != 0)
        return -1;

    p1_im = (log_t *)malloc(len*sizeof(log_t));
    p2_im = (log_t *)malloc(len*sizeof(log_t));
    p3_im = (log_t *)malloc(len*sizeof(log_t));

    for(i = 0; i < len; i++)
    {
        p1_im[i].value = log(0); p1_im[i].sign = +1;
        p2_im[i].value = log(0); p2_im[i].sign = +1;
        p3_im[i].value = log(0); p3_im[i].sign = +1;
    }

    fft(&self, p1, p1_im);
    fft(&self, p2, p2_im);
    fft(&self, p3, p3_im);

    /*
    fft_permute_bitrev(&self, p1, p1_im);
    fft_permute_bitrev(&self, p2, p2_im);
    fft_permute_bitrev(&self, p3, p3_im);
    */

    for(i = 0; i < len; i++)
    {
        double re, im;
        int re_sign, im_sign;
        /* p2*p3 => p2 */
        re = logadd_s(p2[i].value+p3[i].value, p2[i].sign*p3[i].sign, p2_im[i].value+p3_im[i].value, -p2_im[i].sign*p3_im[i].sign, &re_sign);
        im = logadd_s(p2_im[i].value+p3[i].value, p2_im[i].sign*p3[i].sign, p2[i].value+p3_im[i].value, p2[i].sign*p3_im[i].sign, &im_sign);
        p2[i].value    = re;
        p2[i].sign     = re_sign;
        p2_im[i].value = im;
        p2_im[i].sign  = im_sign;

        /* p2*p1 => p1 */
        re = logadd_s(p2[i].value+p1[i].value, p2[i].sign*p1[i].sign, p2_im[i].value+p1_im[i].value, -p2_im[i].sign*p1_im[i].sign, &re_sign);
        im = logadd_s(p2_im[i].value+p1[i].value, p2_im[i].sign*p1[i].sign, p2[i].value+p1_im[i].value, p2[i].sign*p1_im[i].sign, &im_sign);
        p1[i].value    = re;
        p1[i].sign     = re_sign;
        p1_im[i].value = im;
        p1_im[i].sign  = im_sign;

        //printf("%d: M %g%+gi\n", i, p1[i].sign*exp(p1[i].value), p1_im[i].sign*exp(p1_im[i].value));
    }

    fft_permute_bitrev(&self, p1, p1_im);
    ifft(&self, p1, p1_im);
    fft_permute_bitrev(&self, p1, p1_im);
    fft_free(&self);

    free(p1_im);
    free(p2_im);
    free(p3_im);

    return 0;
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
static void _fft(fft_t *self, int dir, log_t *A_re, log_t *A_im) 
{
    const int n = self->n;
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

void ifft(fft_t *self, log_t *A_re, log_t *A_im)
{
    _fft(self, -1, A_re, A_im);

}

void fft(fft_t *self, log_t *A_re, log_t *A_im) 
{
    _fft(self, +1, A_re, A_im);
}

#if 0
int main(int argc, char *argv[])
{
    int n = 4;
    int i;
    log_t p1[n], p2[n], p3[n];
    
    for(i = 0; i < n; i++)
    {
        p1[i].sign = 1;
        p2[i].sign = 1;
        p3[i].sign = 1;

        if(i > 1)
        {
            p1[i].value = log(0);
            p2[i].value = log(0);
            p3[i].value = log(0);
        }
        else
        {
            p1[i].value = log(i+1);
            p2[i].value = log(i+2);
            p3[i].value = log(i+3);
        }
    }

    fft_polymult3(p1, p2, p3, n);
    
    printf("\n\n");
    for(i = 0; i < n; i++)
        printf("%d: %+g\n", i, p1[i].sign*exp(p1[i].value));

    return 0;
}
#endif
