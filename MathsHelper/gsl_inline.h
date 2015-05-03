/* gsl_inline.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.  

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
	#ifdef HAVE_INLINE
	INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif
   
*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */
