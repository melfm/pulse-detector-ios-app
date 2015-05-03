/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

extern double gsl_log1p (const double x);
extern double gsl_expm1 (const double x);
extern double gsl_hypot (const double x, const double y);
extern double gsl_hypot3 (const double x, const double y, const double z);
extern double gsl_acosh (const double x);
extern double gsl_asinh (const double x);
extern double gsl_atanh (const double x);

extern int gsl_isnan (const double x);
extern int gsl_isinf (const double x);
extern int gsl_finite (const double x);

extern double gsl_nan (void);
extern double gsl_posinf (void);
extern double gsl_neginf (void);
extern double gsl_fdiv (const double x, const double y);

extern double gsl_coerce_double (const double x);
extern float gsl_coerce_float (const float x);
extern long double gsl_coerce_long_double (const long double x);

extern double gsl_ldexp(const double x, const int e);
extern double gsl_frexp(const double x, int * e);

extern int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */
