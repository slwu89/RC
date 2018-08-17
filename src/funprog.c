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
  SEXP out  = PROTECT(VECTOR_ELT(x,0));
  Rprintf("out is now: %i\n",asInteger(out));
  pCalls += 1;

  /* make a symbol-value for the "i" part of [[i]] */
  SEXP iSym = PROTECT(install("i"));
  SEXP iVal = PROTECT(ScalarInteger(1));
  defineVar(iSym,iVal,rho);
  pCalls += 2;

  /* make the part to do: f(out, x[[i]]) */
  SEXP bracket = PROTECT(LCONS(R_Bracket2Symbol, LCONS(x, LCONS(iSym, R_NilValue)))); /* x[[i]] */
  SEXP R_fcall = PROTECT(LCONS(f, LCONS(out, LCONS(bracket, R_NilValue))));
  pCalls += 2;

  /* print the crap */
  Rprintf("printing the call we built\n");
  SEXP printSym = PROTECT(install("print"));
  SEXP fcall = PROTECT(LCONS(printSym,LCONS(fcall,R_NilValue)));
  pCalls += 2;
  eval(fcall,rho);

  for(int i=2; i<=n; ++i){
    INTEGER(iVal)[0] = i;
    Rprintf("we are at INTEGER(iVal)[0]: %i\n",INTEGER(iVal)[0]);
    out = eval(R_fcall,rho);
    Rprintf("out is now: %i\n",asInteger(out));
  }

  UNPROTECT(pCalls);
  return out;
};

//out <- x[[1]]
//for(i in seq(2, length(x))) {
//  out <- f(out, x[[i]])
//}
//out
