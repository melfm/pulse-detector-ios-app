/* fft/gsl_fft_halfcomplex.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

#ifndef __GSL_FFT_HALFCOMPLEX_H__
#define __GSL_FFT_HALFCOMPLEX_H__

#include <stddef.h>

#include "gsl_math.h"
#include "gsl_complex.h"
#include "gsl_fft.h"
#include "gsl_fft_real.h"

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

int gsl_fft_halfcomplex_radix2_backward (double data[], const size_t stride, const size_t n);
int gsl_fft_halfcomplex_radix2_inverse (double data[], const size_t stride, const size_t n);
int gsl_fft_halfcomplex_radix2_transform (double data[], const size_t stride, const size_t n);

typedef struct
  {
    size_t n;
    size_t nf;
    size_t factor[64];
    gsl_complex *twiddle[64];
    gsl_complex *trig;
  }
gsl_fft_halfcomplex_wavetable;

gsl_fft_halfcomplex_wavetable * gsl_fft_halfcomplex_wavetable_alloc (size_t n);

void
gsl_fft_halfcomplex_wavetable_free (gsl_fft_halfcomplex_wavetable * wavetable);


int gsl_fft_halfcomplex_backward (double data[], const size_t stride, const size_t n,
                                  const gsl_fft_halfcomplex_wavetable * wavetable,
                                  gsl_fft_real_workspace * work);

int gsl_fft_halfcomplex_inverse (double data[], const size_t stride, const size_t n,
                                 const gsl_fft_halfcomplex_wavetable * wavetable,
                                 gsl_fft_real_workspace * work);

int gsl_fft_halfcomplex_transform (double data[], const size_t stride, const size_t n,
                                   const gsl_fft_halfcomplex_wavetable * wavetable,
                                   gsl_fft_real_workspace * work);

extern int
gsl_fft_halfcomplex_unpack (const double halfcomplex_coefficient[],
                            double complex_coefficient[],
                            const size_t stride, const size_t n);

int
gsl_fft_halfcomplex_radix2_unpack (const double halfcomplex_coefficient[],
                                   double complex_coefficient[],
                                   const size_t stride, const size_t n);

__END_DECLS

#endif /* __GSL_FFT_HALFCOMPLEX_H__ */
