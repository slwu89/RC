/* ##################################################
#   MCR (Mutagenic Chain Reaction) likelihood function
#   Sean Wu
################################################## */

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

/* MCR model likelihood across all experiments */
SEXP C_logLike_AllExpts(SEXP c_R, SEXP pH_R, SEXP prR_R, SEXP sH_R, SEXP sR_R, SEXP sB_R,
                        SEXP GFPp_ym_F_1_R, SEXP GFPp_yp_F_1_R, SEXP GFPm_ym_F_1_R, SEXP GFPm_yp_F_1_R,
                        SEXP GFPp_ym_M_1_R, SEXP GFPp_yp_M_1_R, SEXP GFPm_ym_M_1_R, SEXP GFPm_yp_M_1_R,
                        SEXP GFPp_ym_F_2_R, SEXP GFPp_yp_F_2_R, SEXP GFPm_ym_F_2_R, SEXP GFPm_yp_F_2_R,
                        SEXP GFPp_ym_M_2_R, SEXP GFPp_yp_M_2_R, SEXP GFPm_ym_M_2_R, SEXP GFPm_yp_M_2_R,
                        SEXP GFPp_ym_F_3_R, SEXP GFPp_yp_F_3_R, SEXP GFPm_ym_F_3_R, SEXP GFPm_yp_F_3_R,
                        SEXP GFPp_ym_M_3_R, SEXP GFPp_yp_M_3_R, SEXP GFPm_ym_M_3_R, SEXP GFPm_yp_M_3_R,
                        SEXP GFPp_ym_F_4_R, SEXP GFPp_yp_F_4_R, SEXP GFPm_ym_F_4_R, SEXP GFPm_yp_F_4_R,
                        SEXP GFPp_ym_M_4_R, SEXP GFPp_yp_M_4_R, SEXP GFPm_ym_M_4_R, SEXP GFPm_yp_M_4_R,
                        SEXP gens_R);

#endif
