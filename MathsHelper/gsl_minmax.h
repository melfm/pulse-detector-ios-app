/* gsl_minmax.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

#ifndef __GSL_MINMAX_H__
#define __GSL_MINMAX_H__
#include "gsl_inline.h"

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

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double gsl_max (double a, double b);
double gsl_min (double a, double b);

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

INLINE_FUN int GSL_MAX_INT (int a, int b);
INLINE_FUN int GSL_MIN_INT (int a, int b);
INLINE_FUN double GSL_MAX_DBL (double a, double b);
INLINE_FUN double GSL_MIN_DBL (double a, double b);
INLINE_FUN long double GSL_MAX_LDBL (long double a, long double b);
INLINE_FUN long double GSL_MIN_LDBL (long double a, long double b);

INLINE_FUN int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_POW_INT_H__ */
