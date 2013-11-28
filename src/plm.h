#ifndef __PLM_H
#define __PLM_H

#include <gsl/gsl_matrix.h>

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

double Plm(int l, int m, double x, gsl_matrix **pmPlm, gsl_matrix **pmdPlm);
double dPlm(int l, int m, double x);

double Pl1mPl2m(int l1, int l2, int m, double x);
double dPl1mdPl2m(int l1, int l2, int m, double x);
double Pl1mdPl2m(int l1, int l2, int m, double x);

double Yl1mdYl2m(int l1, int l2, int m, double x);
double dYl1mdYl2m(int l1, int l2, int m, double x);
double Yl1mYl2m(int l1, int l2, int m, double x);
#endif
