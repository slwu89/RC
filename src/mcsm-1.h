/* --------------------------------------------------------------------------------
#   Monte Carlo Statistical Methods: Ch1
 -------------------------------------------------------------------------------- */

#ifndef MCSM_1
#define MCSM_1

/* C headers */
#include <stdio.h>
#include <math.h>

/* R's C API */
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

/* A.1: pool-adjacent violators */
SEXP mcsm_pava_c(SEXP f, SEXP w);

#endif
