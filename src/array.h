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
SEXP sub_array3_C(SEXP iR, SEXP jR, SEXP kR, SEXP arrR);

/* subset a 4d array (4-tensor) */
SEXP sub_array4_C(SEXP iR, SEXP jR, SEXP kR, SEXP lR, SEXP arrR);

#endif
