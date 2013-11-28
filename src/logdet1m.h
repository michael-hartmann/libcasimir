#ifndef __LOGDET1M_H
#define __LOGDET1M_H

int logdet1m_taylor(gsl_matrix *M, double *logdet);
int logdet1m_eigenvalues(gsl_matrix *M, double *logdet);

#endif
