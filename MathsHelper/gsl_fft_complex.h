/* fft/gsl_fft_complex.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

#ifndef __GSL_FFT_COMPLEX_H__
#define __GSL_FFT_COMPLEX_H__

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

/*  Power of 2 routines  */


int gsl_fft_complex_radix2_forward (gsl_complex_packed_array data,
                                    const size_t stride,
                                    const size_t n);

int gsl_fft_complex_radix2_backward (gsl_complex_packed_array data,
                                     const size_t stride,
                                     const size_t n);

int gsl_fft_complex_radix2_inverse (gsl_complex_packed_array data,
                                    const size_t stride,
                                    const size_t n);

int gsl_fft_complex_radix2_transform (gsl_complex_packed_array data,
                                      const size_t stride,
                                      const size_t n,
                                      const gsl_fft_direction sign);

int gsl_fft_complex_radix2_dif_forward (gsl_complex_packed_array data,
                                        const size_t stride,
                                        const size_t n);

int gsl_fft_complex_radix2_dif_backward (gsl_complex_packed_array data,
                                         const size_t stride,
                                         const size_t n);

int gsl_fft_complex_radix2_dif_inverse (gsl_complex_packed_array data,
                                        const size_t stride,
                                        const size_t n);

int gsl_fft_complex_radix2_dif_transform (gsl_complex_packed_array data,
                                          const size_t stride,
                                          const size_t n,
                                          const gsl_fft_direction sign);

/*  Mixed Radix general-N routines  */

typedef struct
  {
    size_t n;
    size_t nf;
    size_t factor[64];
    gsl_complex *twiddle[64];
    gsl_complex *trig;
  }
gsl_fft_complex_wavetable;

typedef struct
{
  size_t n;
  double *scratch;
}
gsl_fft_complex_workspace;


gsl_fft_complex_wavetable *gsl_fft_complex_wavetable_alloc (size_t n);

void gsl_fft_complex_wavetable_free (gsl_fft_complex_wavetable * wavetable);

gsl_fft_complex_workspace *gsl_fft_complex_workspace_alloc (size_t n);

void gsl_fft_complex_workspace_free (gsl_fft_complex_workspace * workspace);

int gsl_fft_complex_memcpy (gsl_fft_complex_wavetable * dest,
                            gsl_fft_complex_wavetable * src);


int gsl_fft_complex_forward (gsl_complex_packed_array data,
                             const size_t stride,
                             const size_t n,
                             const gsl_fft_complex_wavetable * wavetable,
                             gsl_fft_complex_workspace * work);

int gsl_fft_complex_backward (gsl_complex_packed_array data,
                              const size_t stride,
                              const size_t n,
                              const gsl_fft_complex_wavetable * wavetable,
                              gsl_fft_complex_workspace * work);

int gsl_fft_complex_inverse (gsl_complex_packed_array data,
                             const size_t stride,
                             const size_t n,
                             const gsl_fft_complex_wavetable * wavetable,
                             gsl_fft_complex_workspace * work);

int gsl_fft_complex_transform (gsl_complex_packed_array data,
                               const size_t stride, const size_t n,
                               const gsl_fft_complex_wavetable * wavetable,
                               gsl_fft_complex_workspace * work,
                               const gsl_fft_direction sign);

__END_DECLS

#endif /* __GSL_FFT_COMPLEX_H__ */
