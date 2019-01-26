/* ##################################################
#   MGDrivE-Epi (Secret Projects)
#   fixed discrete delay
#   Sean Wu
################################################## */

#include "epi-delay.h"
#include "utilities.h"

/* aquatic population dynamics */
void oneDay_oviposit_C(SEXP cube){

  double time_lag = 0;
  int nG = 0;
  //SEXP timeAq = getListElement(pars,"timeAq");
  //for(size_t i=0; i<length(timeAq); i++){
  //  time_lag += REAL(timeAq)[i];
  //}
  nG = asInteger(getListElement(cube,"genotypesN"));

  printf("time_lag: %d, nG: %d \n",time_lag,nG);
};
