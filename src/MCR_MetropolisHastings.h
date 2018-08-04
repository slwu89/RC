/* ##################################################
#   MCR (Mutagenic Chain Reaction) Metropolis-Hastings MCMC
#   Sean Wu
################################################## */

#ifndef MCR_MCMC_H
#define MCR_MCMC_H

 /* C headers */
#include <stdio.h>
#include <math.h>

 /* R's C API */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* prior function */
SEXP logPrior(SEXP c_R, SEXP pH_R, SEXP prR_R, SEXP sH_R, SEXP sR_R, SEXP sB_R);

#endif
