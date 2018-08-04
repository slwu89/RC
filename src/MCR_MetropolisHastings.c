/* ##################################################
#   MCR (Mutagenic Chain Reaction) Metropolis-Hastings MCMC
#   Sean Wu
################################################## */

#include "MCR_MetropolisHastings.h"
#include "MCR.h"

 /* prior function */
 SEXP C_logPrior(SEXP c_R, SEXP pH_R, SEXP prR_R, SEXP sH_R, SEXP sR_R, SEXP sB_R){
  /* Prior on c: */
  double logPrior_c = dunif(asReal(c_R),0,1,1);
  /* Prior on pH: */
  double logPrior_pH = dunif(asReal(pH_R),0,1,1);
  /* Prior on prR: */
  double logPrior_prR = dunif(asReal(prR_R),0,1,1);
  /* Prior on sH: */
  double logPrior_sH = dunif(asReal(sH_R),-2,1,1);
  /* Prior on sR: */
  double logPrior_sR = dunif(asReal(sR_R),-2,1,1);
  /* Prior on sB: */
  double logPrior_sB = dunif(asReal(sB_R),-2,1,1);

  double out = logPrior_c + logPrior_pH + logPrior_prR + logPrior_prR + logPrior_sH + logPrior_sR + logPrior_sB;
  return ScalarReal(out);
 };


 /* log value of target distribution */
 SEXP C_logPosterior(SEXP c_R, SEXP pH_R, SEXP prR_R, SEXP sH_R, SEXP sR_R, SEXP sB_R,
                   SEXP GFPp_ym_F_1_R, SEXP GFPp_yp_F_1_R, SEXP GFPm_ym_F_1_R, SEXP GFPm_yp_F_1_R,
                   SEXP GFPp_ym_M_1_R, SEXP GFPp_yp_M_1_R, SEXP GFPm_ym_M_1_R, SEXP GFPm_yp_M_1_R,
                   SEXP GFPp_ym_F_2_R, SEXP GFPp_yp_F_2_R, SEXP GFPm_ym_F_2_R, SEXP GFPm_yp_F_2_R,
                   SEXP GFPp_ym_M_2_R, SEXP GFPp_yp_M_2_R, SEXP GFPm_ym_M_2_R, SEXP GFPm_yp_M_2_R,
                   SEXP GFPp_ym_F_3_R, SEXP GFPp_yp_F_3_R, SEXP GFPm_ym_F_3_R, SEXP GFPm_yp_F_3_R,
                   SEXP GFPp_ym_M_3_R, SEXP GFPp_yp_M_3_R, SEXP GFPm_ym_M_3_R, SEXP GFPm_yp_M_3_R,
                   SEXP GFPp_ym_F_4_R, SEXP GFPp_yp_F_4_R, SEXP GFPm_ym_F_4_R, SEXP GFPm_yp_F_4_R,
                   SEXP GFPp_ym_M_4_R, SEXP GFPp_yp_M_4_R, SEXP GFPm_ym_M_4_R, SEXP GFPm_yp_M_4_R,
                   SEXP gens_R){

   /* evaluate prior distribution */
   double logPrior_val = asReal(C_logPrior(c_R,pH_R,prR_R,sH_R,sR_R,sB_R));

   /* evaluate likelihood */
   double logLike_val = asReal(C_logLike_AllExpts(c_R, pH_R, prR_R, sH_R, sR_R, sB_R,
                                                  GFPp_ym_F_1_R, GFPp_yp_F_1_R, GFPm_ym_F_1_R, GFPm_yp_F_1_R,
                                                  GFPp_ym_M_1_R, GFPp_yp_M_1_R, GFPm_ym_M_1_R, GFPm_yp_M_1_R,
                                                  GFPp_ym_F_2_R, GFPp_yp_F_2_R, GFPm_ym_F_2_R, GFPm_yp_F_2_R,
                                                  GFPp_ym_M_2_R, GFPp_yp_M_2_R, GFPm_ym_M_2_R, GFPm_yp_M_2_R,
                                                  GFPp_ym_F_3_R, GFPp_yp_F_3_R, GFPm_ym_F_3_R, GFPm_yp_F_3_R,
                                                  GFPp_ym_M_3_R, GFPp_yp_M_3_R, GFPm_ym_M_3_R, GFPm_yp_M_3_R,
                                                  GFPp_ym_F_4_R, GFPp_yp_F_4_R, GFPm_ym_F_4_R, GFPm_yp_F_4_R,
                                                  GFPp_ym_M_4_R, GFPp_yp_M_4_R, GFPm_ym_M_4_R, GFPm_yp_M_4_R,
                                                  gens_R));
   return ScalarReal(logPrior_val + logLike_val);
 };
