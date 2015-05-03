/* fft/gsl_fft_real.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

#ifndef __GSL_FFT_REAL_H__
#define __GSL_FFT_REAL_H__

#include <stddef.h>

#include "gsl_math.h"
#include "gsl_complex.h"
#include "gsl_fft.h"

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

int gsl_fft_real_radix2_transform (double data[], const size_t stride, const size_t n) ;

typedef struct
  {
    size_t n;
    size_t nf;
    size_t factor[64];
    gsl_complex *twiddle[64];
    gsl_complex *trig;
  }
gsl_fft_real_wavetable;

typedef struct
  {
    size_t n;
    double *scratch;
  }
gsl_fft_real_workspace;

extern gsl_fft_real_wavetable * gsl_fft_real_wavetable_alloc (size_t n);

extern void  gsl_fft_real_wavetable_free (gsl_fft_real_wavetable * wavetable);

extern gsl_fft_real_workspace * gsl_fft_real_workspace_alloc (size_t n);

extern void  gsl_fft_real_workspace_free (gsl_fft_real_workspace * workspace);


extern int gsl_fft_real_transform (double data[], const size_t stride, const size_t n,
                            const gsl_fft_real_wavetable * wavetable,
                            gsl_fft_real_workspace * work);


int gsl_fft_real_unpack (const double real_coefficient[],
                         double complex_coefficient[],
                         const size_t stride, const size_t n);

__END_DECLS

#endif /* __GSL_FFT_REAL_H__ */
