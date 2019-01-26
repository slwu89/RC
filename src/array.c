/* ##################################################
#   interacting with arrays tidbits
#   Sean Wu
################################################## */

#include "array.h"

/* subset a 3d array (3-tensor) */
SEXP sub_array3_C(SEXP iR, SEXP jR, SEXP kR, SEXP arrR){
  SEXP dims;
  int i,j,k;
  i = asInteger(iR)-1;
  j = asInteger(jR)-1;
  k = asInteger(kR)-1;

  dims = getAttrib(arrR, R_DimSymbol);
  int nx = INTEGER(dims)[0];
  int ny = INTEGER(dims)[1];
  int nz = INTEGER(dims)[2];

  if(i < 0 || i >= nx){
    Rprintf("bad value for i index\n");
    return R_NilValue;
  }
  if(j < 0 || j >= ny){
    Rprintf("bad value for j index\n");
    return R_NilValue;
  }
  if(k < 0 || k >= nz){
    Rprintf("bad value for k index\n");
    return R_NilValue;
  }

  int *arr = INTEGER(arrR);
  int ans;
  ans = arr[i + (j * nx) + (k * nx * ny)];
  return ScalarInteger(ans);
};


/* subset a 4d array (4-tensor) */
SEXP sub_array4_C(SEXP iR, SEXP jR, SEXP kR, SEXP lR, SEXP arrR){

  SEXP dims;
  int i,j,k,l;
  i = asInteger(iR)-1;
  j = asInteger(jR)-1;
  k = asInteger(kR)-1;
  l = asInteger(lR)-1;

  dims = getAttrib(arrR, R_DimSymbol);
  int nx = INTEGER(dims)[0];
  int ny = INTEGER(dims)[1];
  int nz0 = INTEGER(dims)[2];
  int nz1 = INTEGER(dims)[4];

  if(i < 0 || i >= nx){
    Rprintf("bad value for i index\n");
    return R_NilValue;
  }
  if(j < 0 || j >= ny){
    Rprintf("bad value for j index\n");
    return R_NilValue;
  }
  if(k < 0 || k >= nz0){
    Rprintf("bad value for k index\n");
    return R_NilValue;
  }
  if(l < 0 || l >= nz1){
    Rprintf("bad value for l index\n");
    return R_NilValue;
  }

  int *arr = INTEGER(arrR);
  int ans;
  ans = arr[i + (j * nx) + (k * nx * ny) + (l * nx * ny * nz0)];
  return ScalarInteger(ans);
};
