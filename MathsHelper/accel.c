/* interpolation/accel.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

/* Author:  G. Jungman
 */
#include "config.h"
#include <stdlib.h>
#include "gsl_errno.h"
#include "gsl_interp.h"

gsl_interp_accel *
gsl_interp_accel_alloc (void)
{
  gsl_interp_accel *a = (gsl_interp_accel *) malloc (sizeof (gsl_interp_accel));
  if (a == 0)
    {
      GSL_ERROR_NULL("could not allocate space for gsl_interp_accel", GSL_ENOMEM);
    }

  a->cache = 0;
  a->hit_count = 0;
  a->miss_count = 0;

  return a;
}

int
gsl_interp_accel_reset (gsl_interp_accel * a)
{
  a->cache = 0;
  a->hit_count = 0;
  a->miss_count = 0;

  return GSL_SUCCESS;
}

void
gsl_interp_accel_free (gsl_interp_accel * a)
{
  RETURN_IF_NULL (a);
  free (a);
}
