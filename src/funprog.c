/* ##################################################
#   Functional Programming tidbits
#   Sean Wu
################################################## */

#include "funprog.h"

/* simple version of Reduce (see http://adv-r.had.co.nz/Functionals.html#functionals-fp) */
SEXP Reduce_Simple_C(SEXP f, SEXP x, SEXP rho){

  size_t pCalls = 0;

  /* allocate output vector */
  size_t n = length(x);

  /* this part: out <- x[[1]] */
  SEXP out =  PROTECT(VECTOR_ELT(x,0));
  pCalls += 1;

  /* make a symbol-value for the "i" part of [[i]] */
  SEXP ii = PROTECT(ScalarInteger(1));
  pCalls += 1;

  /* make the part to do: f(out, x[[i]]) */
  SEXP bracket = PROTECT(LCONS(R_Bracket2Symbol, LCONS(x, LCONS(ii, R_NilValue)))); /* x[[i]] */
  SEXP R_fcall = PROTECT(LCONS(f, LCONS( out, LCONS(bracket, R_NilValue))));
  pCalls += 2;

  /* reduce the list */
  for(int i=2; i<=n; ++i){
    INTEGER(ii)[0] = i;
    out = eval(R_fcall,rho);
    SETCADR(R_fcall,out);
  }

  UNPROTECT(pCalls);
  return out;
};
