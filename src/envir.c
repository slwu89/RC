/* ##################################################
#   interacting with environments tidbits
#   Sean Wu
################################################## */

#include "envir.h"

/* apply a function f to all elements in a hashed envir for side-effects (no return value) */
SEXP eapply_fast_C(SEXP e, SEXP f, SEXP rho){

  int pCalls = 0;

  /* env is the hash table */
  SEXP env = PROTECT(eval(e,rho));
  int n = HashTableSize(HASHTAB(env), 0);
  pCalls += 1;

  /* get the values out of the hash table as vector/list: vals */
  SEXP vals;
  PROTECT(vals = allocVector(VECSXP,n));
  int ix = 0;
  HashTableValues(HASHTAB(env), 0, vals, &ix);
  pCalls += 1;

  /* make the bit of the function call that indexes over the values in the hash table */
  SEXP i = ScalarInteger(1);
  int* iptr = INTEGER(i);
  SEXP bracket;
  PROTECT(bracket = LCONS(R_Bracket2Symbol, LCONS(vals, LCONS(i, R_NilValue))));
  pCalls += 1;

  /* make the f(X[[i]],...) bit of the function call (tmp is the indexing) */
  SEXP R_fcall;
  PROTECT(R_fcall = LCONS(f, LCONS(bracket, LCONS(R_DotsSymbol, R_NilValue))));
  pCalls += 1;

  /* map the function(...) over the hash table (note using R's 1-based indexing) */
  for(int j=1; j<=n; j++){
    *iptr = j;
    eval(R_fcall, rho);
  };

  UNPROTECT(pCalls);
  return R_NilValue;
};
