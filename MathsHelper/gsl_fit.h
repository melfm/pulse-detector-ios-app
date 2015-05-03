/* fit/gsl_fit.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

#ifndef __GSL_FIT_H__
#define __GSL_FIT_H__

#include <stdlib.h>
#include "gsl_math.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

int gsl_fit_linear (const double * x, const size_t xstride,
                    const double * y, const size_t ystride,
                    const size_t n,
                    double * c0, double * c1, 
                    double * cov00, double * cov01, double * cov11, 
                    double * sumsq);


int gsl_fit_wlinear (const double * x, const size_t xstride,
                     const double * w, const size_t wstride,
                     const double * y, const size_t ystride,
                     const size_t n,
                     double * c0, double * c1, 
                     double * cov00, double * cov01, double * cov11, 
                     double * chisq);

int
gsl_fit_linear_est (const double x, 
                    const double c0, const double c1, 
                    const double cov00, const double cov01, const double cov11,
                    double *y, double *y_err);


int gsl_fit_mul (const double * x, const size_t xstride,
                 const double * y, const size_t ystride,
                 const size_t n,
                 double * c1, 
                 double * cov11, 
                 double * sumsq);

int gsl_fit_wmul (const double * x, const size_t xstride,
                  const double * w, const size_t wstride,
                  const double * y, const size_t ystride,
                  const size_t n,
                  double * c1, 
                  double * cov11, 
                  double * sumsq);


int
gsl_fit_mul_est (const double x, 
                 const double c1, 
                 const double cov11,
                 double *y, double *y_err);

__END_DECLS

#endif /* __GSL_FIT_H__ */
