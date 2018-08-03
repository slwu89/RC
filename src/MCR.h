#ifndef MCR_H
#define MCR_H

/* C headers */
#include <stdio.h>
#include <math.h>

/* R's C API */
#include <R.h>
#include <Rinternals.h>

/* MCR model likelihood */
SEXP C_logLike_MCRMod(SEXP c_R, SEXP pH_R, SEXP prR_R, SEXP sH_R, SEXP sR_R, SEXP sB_R,
                      SEXP GFPp_ym_F_R, SEXP GFPp_yp_F_R, SEXP GFPm_ym_F_R, SEXP GFPm_yp_F_R,
                      SEXP GFPp_ym_M_R, SEXP GFPp_yp_M_R, SEXP GFPm_ym_M_R, SEXP GFPm_yp_M_R,
                      SEXP gens_R);

#endif
