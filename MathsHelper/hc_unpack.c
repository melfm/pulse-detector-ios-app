/* fft/hc_unpack.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

#include "gsl_complex.h"
#include "gsl_errno.h"
#include "complex_internal.h"

#include "gsl_fft_halfcomplex.h"

#define ATOMIC double

int gsl_fft_halfcomplex_unpack (const double halfcomplex_coefficient[],
                                      double complex_coefficient[],
                                      const size_t stride, const size_t n)
{
  size_t i;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  REAL(complex_coefficient,stride,0) = halfcomplex_coefficient[0];
  IMAG(complex_coefficient,stride,0) = 0.0;

  for (i = 1; i < n - i; i++)
    {
      const ATOMIC hc_real = halfcomplex_coefficient[(2 * i - 1) * stride];
      const ATOMIC hc_imag = halfcomplex_coefficient[2 * i * stride];

      REAL(complex_coefficient,stride,i) = hc_real;
      IMAG(complex_coefficient,stride,i) = hc_imag;
      REAL(complex_coefficient,stride,n - i) = hc_real;
      IMAG(complex_coefficient,stride,n - i) = -hc_imag;
    }

  if (i == n - i)
    {
      REAL(complex_coefficient,stride,i) = halfcomplex_coefficient[(n - 1) * stride];
      IMAG(complex_coefficient,stride,i) = 0.0;
    }

  return 0;
}


int gsl_fft_halfcomplex_radix2_unpack (const double halfcomplex_coefficient[],
                                             double complex_coefficient[],
                                             const size_t stride, const size_t n)
{
  size_t i;

  if (n == 0)
    {
      GSL_ERROR ("length n must be positive integer", GSL_EDOM);
    }

  REAL(complex_coefficient,stride,0) = halfcomplex_coefficient[0];
  IMAG(complex_coefficient,stride,0) = 0.0;

  for (i = 1; i < n - i; i++)
    {
      const ATOMIC hc_real = halfcomplex_coefficient[i * stride];
      const ATOMIC hc_imag = halfcomplex_coefficient[(n - i) * stride];

      REAL(complex_coefficient,stride,i) = hc_real;
      IMAG(complex_coefficient,stride,i) = hc_imag;
      REAL(complex_coefficient,stride,n - i) = hc_real;
      IMAG(complex_coefficient,stride,n - i) = -hc_imag;
    }

  if (i == n - i)
    {
      REAL(complex_coefficient,stride,i) = halfcomplex_coefficient[i * stride];
      IMAG(complex_coefficient,stride,i) = 0.0;
    }

  return 0;
}

