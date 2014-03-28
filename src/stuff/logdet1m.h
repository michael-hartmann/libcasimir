#ifndef __LOGDET1M_H
#define __LOGDET1M_H

#include <gsl/gsl_matrix.h>

int logdet1m_eigenvalues(gsl_matrix *M, double *logdet);

#endif
