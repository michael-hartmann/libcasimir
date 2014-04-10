#ifndef __POLY_H
#define __POLY_H

void poly_multiply_naive(edouble p1[], size_t len_p1, edouble p2[], size_t len_p2, edouble pdest[]);
void poly_add(edouble p1[], size_t len_p1, edouble p2[], size_t len_p2, edouble pdest[]);
void poly_mult_karatsuba(edouble p1[], size_t len_p1, edouble p2[], size_t len_p2, edouble pdest[]);

#endif
