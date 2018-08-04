/* ##################################################
#   MCR (Mutagenic Chain Reaction) Metropolis-Hastings MCMC
#   Sean Wu
################################################## */

#include "MCR_MetropolisHastings.h"

 /* prior function */
 SEXP logPrior(SEXP c_R, SEXP pH_R, SEXP prR_R, SEXP sH_R, SEXP sR_R, SEXP sB_R){
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
