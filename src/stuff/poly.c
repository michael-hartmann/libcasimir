#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "edouble.h"
#include "sfunc.h"
#include "poly.h"


inline void poly_multiply_naive(edouble p1[], size_t len_p1, edouble p2[], size_t len_p2, edouble pdest[])
{
    size_t i,j;

    for(i = 0; i < len_p1; i++)
        for(j = 0; j < len_p2; j++)
            pdest[i+j] += p1[i]*p2[j];
}

inline void poly_add(edouble p1[], size_t len_p1, edouble p2[], size_t len_p2, edouble pdest[])
{
    size_t i, min = MIN(len_p1, len_p2);

    for(i = 0; i < min; i++)
        pdest[i] = p1[i]+p2[i];
    for(i = min; i < len_p1; i++)
        pdest[i] = p1[i];
    for(i = min; i < len_p2; i++)
        pdest[i] = p2[i];
}

void poly_mult_karatsuba(edouble p1[], size_t len_p1, edouble p2[], size_t len_p2, edouble pdest[])
{
    size_t m2, i, len, len_low1, len_low2, len_high1, len_high2, len_z0, len_z1, len_z2, len_sum1, len_sum2;
    edouble *low1, *low2, *high1, *high2;

    for(i = 0; i < len_p1+len_p2-1; i++)
        pdest[i] = 0;

    if(len_p1 == 1 || len_p2 == 1 || len_p1*len_p2 < 400)
    {
        poly_multiply_naive(p1, len_p1, p2, len_p2, pdest);
        return;
    }

    /* calculates the size of the numbers */
    m2 = MIN(len_p1, len_p2)/2;

    /* split the digit sequences about the middle */
    len_low1  = m2;
    len_high1 = len_p1-m2;
    low1  = p1;
    high1 = p1+m2;

    len_low2  = m2;
    len_high2 = len_p2-m2;
    low2  = p2;
    high2 = p2+m2;

    len_z0 = len_low1+len_low2-1;
    len_z1 = MAX(len_low1,len_low2)+MAX(len_low2, len_high2)-1;
    len_z2 = len_high1+len_high2-1;

    len_sum1 = MAX(len_low1, len_high1);
    len_sum2 = MAX(len_low2, len_high2);

    len = MAX(MAX(len_z0, len_z1), len_z2);
    edouble z0[len];
    edouble z1[len];
    edouble z2[len];
    edouble sum1[len_sum1];
    edouble sum2[len_sum2];

    for(i = 0; i < len; i++)
        z0[i] = z1[i] = z2[i] = 0;

    poly_add(low1, len_low1, high1, len_high1, sum1);
    poly_add(low2, len_low2, high2, len_high2, sum2);

    /* 3 calls made to numbers approximately half the size */
    poly_mult_karatsuba(low1,  len_low1,  low2,  len_low2,  z0);
    poly_mult_karatsuba(sum1,  len_sum1,  sum2,  len_sum2,  z1);
    poly_mult_karatsuba(high1, len_high1, high2, len_high2, z2);

    for(i = 0; i < len; i++)
    {
        pdest[i]      += z0[i];
        pdest[i+m2]   += z1[i]-z2[i]-z0[i];
        pdest[i+2*m2] += z2[i];
    }

    return;
}

#if 0
int main(int argc, char *argv[])
{
    size_t i;
    double p1[] = {1,2,7,19,0,4,-1,-1};
    double p2[] = {3,-4,5,-3};
    double p3[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    poly_mult_karatsuba(p1,8,p2,4,p3);

    printf("\n");
    for(i = 0; i < 11; i++)
        printf("x^%d: %g\n", (int)i, p3[i]);

    return 0;
}
#endif
