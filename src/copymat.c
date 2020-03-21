/* ##################################################
#   testing how to copy subsets of matrices
#   Sean Wu
################################################## */

#include "copymat.h"

/* subset matrix and return */
SEXP submat_C(SEXP mat, SEXP nrow, SEXP ncols){

  int n = Rf_asInteger(nrow);
  int y = Rf_asInteger(ncols);

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP,n,y));

  memmove(REAL(out),REAL(mat),n*y*sizeof(double));

  UNPROTECT(1);
  return out;
};


SEXP submat_cols_C(SEXP mat, SEXP nrows, SEXP ncol){

  int n = Rf_asInteger(nrows);
  int y = Rf_asInteger(ncol);

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP,n,y));

  memmove(REAL(out),REAL(mat),n*y*sizeof(double));

  UNPROTECT(1);
  return out;
};
