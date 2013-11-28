#ifndef __GAUSSLAGUERRE_H
#define __GAUSSLAGUERRE_H


double integrate_gausslaguerre(double(callback(double,void*)), void *params, int n);
void integrateABCD_gausslaguerre(void(callback(double,void*,double *,double *, double *, double *)), void *params, int n, double *A, double *B, double *C, double *D);

#endif
