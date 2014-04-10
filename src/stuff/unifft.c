#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>   

#include "edouble.h"
#include "unifft.h"


/* treats inp as a numbits number and bitreverses it. 
 * inp < 2^(numbits) for meaningful bit-reversal
 */ 
static int bitrev(int inp, int numbits)
{
    int i, rev = 0;
    for(i = 0; i < numbits; i++)
    {
        rev = (rev << 1) | (inp & 1);
        inp >>= 1;
    }
    return rev;
}

int fft_polymult(edouble p1[], edouble p2[], size_t len)
{
    int i;
    fft_t self;
    edouble *p1_im = NULL, *p2_im = NULL;

    p1_im = (edouble *)malloc(len*sizeof(edouble));
    p2_im = (edouble *)malloc(len*sizeof(edouble));

    if(p1_im == NULL || p2_im == NULL)
    {
        if(p1_im != NULL)
            free(p1_im);
        if(p2_im != NULL)
            free(p2_im);

        return 0;
    }

    for(i = 0; i < len; i++)
        p1_im[i] = p2_im[i] = 0;

    if(fft_init(&self, len) == 0)
        return 0;

    fft_fft(&self, 1, p1, p1_im);
    fft_fft(&self, 1, p2, p2_im);

    for(i = 0; i < len; i++)
    {
        edouble re = p1[i]*p2[i]-p1_im[i]*p2_im[i];
        edouble im = p1[i]*p2_im[i]+p1_im[i]*p2[i];

        p1[i]    = re;
        p1_im[i] = im;
    }

    fft_permute_bitrev(&self, p1, p1_im);

    fft_fft(&self, -1, p1, p1_im);
    fft_permute_bitrev(&self, p1, p1_im);

    fft_free(&self);

    return 1;
}

/* returns log n (to the base 2), if n is positive and power of 2 */ 
static int log_2(int n) 
{
    int res; 
    for (res=0; n >= 2; res++) 
        n = n >> 1; 
    return res; 
}
 

/* W will contain roots of unity so that W[bitrev(i,log2n-1)] = e^(2*pi*i/n)
 * n should be a power of 2
 * Note: W is bit-reversal permuted because fft(..) goes faster if this is done.
 *       see that function for more details on why we treat 'i' as a (log2n-1) bit number.
 */
int fft_init(fft_t *self, int n)
{
    int i, br;
    int log2n = log_2(n);

    self->n     = n;
    self->log2n = log2n;

    self->W_re = (edouble *)malloc((n/2)*sizeof(edouble));
    self->W_im = (edouble *)malloc((n/2)*sizeof(edouble));

    if(self->W_re == NULL || self->W_im == NULL)
    {
        fft_free(self);
        return 0;
    }

    for (i=0; i<(n/2); i++)
    {
        br = bitrev(i,log2n-1); 
        self->W_re[br] = +cosq(i*2.0*M_PI/n);
        self->W_im[br] = -sinq(i*2.0*M_PI/n);
    }

    return 1;
}

void fft_free(fft_t *self)
{
    if(self->W_re != NULL)
        free(self->W_re);
    if(self->W_im != NULL)
        free(self->W_im);

    self->W_im = self->W_re = NULL; 
}

#if 0
/* gets no. of points from the user, initialize the points and roots of unity lookup table 
 * and lets fft go. finally bit-reverses the results and outputs them into a file. 
 * n should be a power of 2. 
 */ 
int main(int argc, char *argv[])
{
  int i,n;
  fft_t self;
  edouble A_re[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  edouble A_im[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
  edouble x[] = {1, 2, -4, 0, 0, 0, 0, 0 };
  edouble y[] = {-5, 3, -1, 0, 0, 0, 0, 0 };
  
  n = 8;

  fft_init(&self, n);
  
  fft_fft(&self, 1, A_re, A_im);
  fft_permute_bitrev(&self, A_re, A_im);

  fft_free(&self);

  for(i = 0; i < n; i++)
      printf("%g%+g\n", (double)A_re[i], (double)A_im[i]);
    
  printf("\n");
  fft_polymult(x,y,8);
  for(i = 0; i < n; i++)
      printf("%f*x^%d\n", (double)x[i], i);
  return 0;
}
#endif


/* permutes the array using a bit-reversal permutation */ 
void fft_permute_bitrev(fft_t *self, edouble *A_re, edouble *A_im) 
{ 
  int i, bri;
  edouble t_re, t_im;

  int n     = self->n;
  int log2n = self->log2n;
  
  for (i=0; i<n; i++)
  {
      bri = bitrev(i, log2n);

      /* skip already swapped elements */
      if (bri <= i) continue;

      t_re = A_re[i];
      t_im = A_im[i];
      A_re[i]= A_re[bri];
      A_im[i]= A_im[bri];
      A_re[bri]= t_re;
      A_im[bri]= t_im;
  }  
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
int fft_fft(fft_t *self, int dir, edouble *A_re, edouble *A_im) 
{
    edouble w_re, w_im, u_re, u_im, t_re, t_im;
    int m, g, b;
    int mt, k;
    int n = self->n;
    edouble *W_re = self->W_re;
    edouble *W_im = self->W_im;

    if(W_re == NULL || W_im == NULL)
        return 0;
    
    /* for each stage */  
    for(m=n; m>=2; m=m>>1) 
    {
        /* m = n/2^s; mt = m/2; */
        mt = m >> 1;
        
        /* for each group of butterfly */ 
        for (g=0,k=0; g<n; g+=m,k++) 
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

    if(dir < 0)
    {
        for(m = 0; m < n; m++)
        {
            A_re[m] /= n;
            A_im[m] /= n;
        }
    }

    return 1;
}
