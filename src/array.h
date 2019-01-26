/* ##################################################
#   interacting with arrays tidbits
#   Sean Wu
################################################## */

#ifndef ENVIR_H
#define ENVIR_H

/* R's C API headers */
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

/* C headers */
#include <stdio.h>

/* subset a 3d array (3-tensor) */
SEXP sub_array3_C(SEXP i, SEXP j, SEXP k, SEXP arr);


#endif
