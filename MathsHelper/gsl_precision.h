/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
#include "gsl_types.h"

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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */
