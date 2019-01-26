/* ##################################################
#   Utilities
#   Sean Wu
################################################## */

#include "utilities.h"

/* https://github.com/sbfnk/dynmod/blob/master/src/utilities.c
 * get an element of a list, return null if it doesn't exist
 */
SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  for (i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }

    if (elmt == R_NilValue)
      return NULL;
    else
      return elmt;
}

/* get named element of environment */
SEXP getvar(SEXP name, SEXP rho){
  SEXP ans;

  if(!isString(name) || length(name) != 1)
    error("name is not a single string");
  if(!isEnvironment(rho))
    error("rho should be an environment");
  ans = findVar(installChar(STRING_ELT(name, 0)), rho);
  return ans;
}
