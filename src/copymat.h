/* ##################################################
#   testing how to copy subsets of matrices
#   Sean Wu
################################################## */

#ifndef COPYMAT_H
#define COPYMAT_H

/* R's C API headers */
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

/* C headers */
#include <stdio.h>
#include <stdlib.h>

/* subset matrix and return */
SEXP submat_C(SEXP mat, SEXP nrow, SEXP ncols);

SEXP submat_cols_C(SEXP mat, SEXP nrows, SEXP ncol);


#endif
