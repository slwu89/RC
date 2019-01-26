/* ##################################################
#   Utilities
#   Sean Wu
################################################## */

#ifndef UTILITIES_H
#define UTILITIES_H

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <math.h>

/* https://github.com/sbfnk/dynmod/blob/master/src/utilities.c */
SEXP getListElement(SEXP list, const char *str);

#endif
