#ifndef __GAUSSLAGUERRE_H
#define __GAUSSLAGUERRE_H

double gausslaguerre_integrate(double(callback(double,void*)), void *params, int n);
void   gausslaguerre_integrate_vec(void(callback(double,void*,double *,int)), void *params, int n, double *vec, int len);

#endif
