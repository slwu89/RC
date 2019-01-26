/* ##################################################
#   MGDrivE-Epi (Secret Projects)
#   fixed discrete delay
#   Sean Wu
################################################## */

#ifndef EPI_DELAY_H
#define EPI_DELAY_H

 /* C headers */
#include <stdio.h>

 /* R's C API */
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

/* aquatic population dynamics */
void oneDay_oviposit_C(SEXP cube);

#endif
