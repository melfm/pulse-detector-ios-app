/* err/error.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 *
 */

#include "config.h"
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include "gsl_errno.h"
#include "gsl_message.h"

gsl_error_handler_t * gsl_error_handler = NULL;

static void no_error_handler (const char *reason, const char *file, int line, int gsl_errno);

void
gsl_error (const char * reason, const char * file, int line, int gsl_errno)
{
  if (gsl_error_handler) 
    {
      (*gsl_error_handler) (reason, file, line, gsl_errno);
      return ;
    }

  gsl_stream_printf ("ERROR", file, line, reason);

  fflush (stdout);
  fprintf (stderr, "Default GSL error handler invoked.\n");
  fflush (stderr);

  abort ();
}

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler)
{
  gsl_error_handler_t * previous_handler = gsl_error_handler;
  gsl_error_handler = new_handler;
  return previous_handler;
}


gsl_error_handler_t *
gsl_set_error_handler_off (void)
{
  gsl_error_handler_t * previous_handler = gsl_error_handler;
  gsl_error_handler = no_error_handler;
  return previous_handler;
}

static void
no_error_handler (const char *reason, const char *file, int line, int gsl_errno)
{
  /* do nothing */
  reason = 0;
  file = 0;
  line = 0;
  gsl_errno = 0;
  return;
}


