/* --------------------------------------------------------------------------------
#   Monte Carlo Statistical Methods: Ch1
 -------------------------------------------------------------------------------- */

#include "mcsm-1.h"

/* A.1: pool-adjacent violators */
SEXP mcsm_pava_c(SEXP f, SEXP w){

  int n = length(f);
  if(n != length(w)){
    error(" --- f and w must be equal length --- \n");
  }
  return R_NilValue;
};
