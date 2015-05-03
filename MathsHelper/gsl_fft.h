/* fft/gsl_fft.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

#ifndef __GSL_FFT_H__
#define __GSL_FFT_H__

#include <gsl/gsl_complex.h>

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

typedef enum
  {
    gsl_fft_forward = -1, gsl_fft_backward = +1   
  }
gsl_fft_direction;

/* this gives the sign in the formula

   h(f) = \sum x(t) exp(+/- 2 pi i f t) 
       
   where - is the forward transform direction and + the inverse direction */

__END_DECLS

#endif /* __GSL_FFT_H__ */
