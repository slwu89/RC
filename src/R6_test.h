/* ##################################################
#   Testing C code for R6 objects
#   Sean Wu
################################################## */

#ifndef R6_H
#define R6_H

 /* C headers */
#include <stdio.h>

 /* R's C API */
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

/* simple function to find things in R6 objects */
SEXP find_var_C(SEXP env, SEXP sym, SEXP rho);

#endif
