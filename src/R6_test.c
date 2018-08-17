/* ##################################################
#   Testing C code for R6 objects
#   Sean Wu
################################################## */

#include "R6_test.h"

 /* simple function to find public things in environments and print them */
SEXP find_var_C(SEXP env, SEXP sym, SEXP rho){

  size_t pCalls = 0;

  SEXP val = PROTECT(findVar(installChar(STRING_ELT(sym,0)),env));
  pCalls += 1;

  /* print it */
  SEXP fn = PROTECT(install("print"));
  SEXP fcall = PROTECT(LCONS(fn,LCONS(val,R_NilValue)));
  pCalls += 2;
  eval(fcall,rho);

  UNPROTECT(pCalls);
  return R_NilValue;
};
