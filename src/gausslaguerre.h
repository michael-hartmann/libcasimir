#ifndef __GAUSSLAGUERRE_H
#define __GAUSSLAGUERRE_H

double gauss_laguerre_integrate(double(f(double,void*)), void *params);
void   gausslaguerre_integrate_vec(void(callback(double,void*,double *,int)), void *params, double *vec, int len);

#endif
