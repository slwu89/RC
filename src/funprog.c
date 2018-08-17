/* ##################################################
#   Functional Programming tidbits
#   Sean Wu
################################################## */

#include "funprog.h"

/* simple version of Reduce (see http://adv-r.had.co.nz/Functionals.html#functionals-fp) */
SEXP Reduce_Simple_C(SEXP f, SEXP x, SEXP rho){

  size_t pCalls = 0;

  /* index of where we are in the input vector(list) */
  //SEXP ind = PROTECT(allocVector(INTSXP, 1));
  //pCalls += 1;

  /* allocate output vector */
  size_t n = length(x);

  /* this part: out <- x[[1]] */
  SEXP outSym  = PROTECT(install("out"));
  SEXP outVal =  PROTECT(VECTOR_ELT(x,0));
  defineVar(outSym,outVal,rho);
  Rprintf("out is now: %i\n",asInteger(outVal));
  pCalls += 2;

  /* make a symbol-value for the "i" part of [[i]] */
  SEXP iSym = PROTECT(install("i"));
  SEXP iVal = PROTECT(ScalarInteger(1));
  defineVar(iSym,iVal,rho);
  pCalls += 2;

  /* make the part to do: f(out, x[[i]]) */
  SEXP bracket = PROTECT(LCONS(R_Bracket2Symbol, LCONS(x, LCONS(iSym, R_NilValue)))); /* x[[i]] */
  SEXP R_fcall = PROTECT(LCONS(f, LCONS( outSym, LCONS(bracket, R_NilValue))));
  pCalls += 2;

  Rprintf("starting the loop\n");
  for(int i=1; i<=n; ++i){
    INTEGER(iVal)[0] = i;
    defineVar(iSym,iVal,rho);
    Rprintf("we are at INTEGER(iVal)[0]: %i\n",INTEGER(iVal)[0]);
    outVal = eval(R_fcall,rho);
    Rprintf("out is now: %i\n",asInteger(outVal));
  }

  UNPROTECT(pCalls);
  return outVal;
};

//out <- x[[1]]
//for(i in seq(2, length(x))) {
//  out <- f(out, x[[i]])
//}
//out
