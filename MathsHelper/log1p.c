/* sys/log1p.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

#include "config.h"
#include <math.h>
#include "gsl_sys.h"

double gsl_log1p (const double x)
{
  volatile double y, z;
  y = 1 + x;
  z = y - 1;
  return log(y) - (z-x)/y ;  /* cancels errors with IEEE arithmetic */
}
