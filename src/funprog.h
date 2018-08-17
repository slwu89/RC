/* ##################################################
#   Functional Programming tidbits
#   Sean Wu
################################################## */

#ifndef FUNPROG_H
#define FUNPROG_H

 /* C headers */
#include <stdio.h>

 /* R's C API */
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

/* simple version of Reduce (see http://adv-r.had.co.nz/Functionals.html#functionals-fp) */
SEXP Reduce_Simple_C(SEXP f, SEXP x, SEXP rho);

#endif
