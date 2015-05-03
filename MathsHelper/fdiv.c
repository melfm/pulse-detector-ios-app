/* sys/fdiv.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

#include "config.h"
#include <math.h>
#include "gsl_sys.h"

double 
gsl_fdiv (const double x, const double y)
{
  return x / y;
}
