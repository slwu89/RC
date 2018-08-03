#include "MCR.h"

/* MCR model likelihood */
SEXP C_logLike_MCRMod(SEXP c_R, SEXP pH_R, SEXP prR_R, SEXP sH_R, SEXP sR_R, SEXP sB_R,
                      SEXP GFPp_ym_F_R, SEXP GFPp_yp_F_R, SEXP GFPm_ym_F_R, SEXP GFPm_yp_F_R,
                      SEXP GFPp_ym_M_R, SEXP GFPp_yp_M_R, SEXP GFPm_ym_M_R, SEXP GFPm_yp_M_R,
                      SEXP gens_R)
{

  /* model parameters as C types */
  double c = asReal(c_R);
  double pH = asReal(pH_R);
  double prR = asReal(prR_R);
  double sH = asReal(sH_R);
  double sR = asReal(sR_R);
  double sB = asReal(sB_R);

  int gens = asInteger(gens_R);

  /* derived parameters */
  double pR = (1.0 - pH)*prR;
  double pB = (1.0 - pH)*(1.0 - prR);

  /* vectors of total populations */
  int flen = length(GFPp_ym_F_R);
  int Total_F[flen];
  int* GFPp_ym_F_R_ptr = INTEGER(GFPp_ym_F_R);
  int* GFPp_yp_F_R_ptr = INTEGER(GFPp_yp_F_R);
  int* GFPm_ym_F_R_ptr = INTEGER(GFPm_ym_F_R);
  int* GFPm_yp_F_R_ptr = INTEGER(GFPm_yp_F_R);
  for(int i=0; i <flen; i++){
    Total_F[i] = GFPp_ym_F_R_ptr[i] + GFPp_yp_F_R_ptr[i] + GFPm_ym_F_R_ptr[i] + GFPm_yp_F_R_ptr[i];
  }

  int mlen = length(GFPp_ym_M_R);
  int Total_M[mlen];
  int* GFPp_ym_M_R_ptr = INTEGER(GFPp_ym_M_R);
  int* GFPp_yp_M_R_ptr = INTEGER(GFPp_yp_M_R);
  int* GFPm_ym_M_R_ptr = INTEGER(GFPm_ym_M_R);
  int* GFPm_yp_M_R_ptr = INTEGER(GFPm_yp_M_R);
  for(int i=0; i <mlen; i++){
    Total_M[i] = GFPp_ym_M_R_ptr[i] + GFPp_yp_M_R_ptr[i] + GFPm_ym_M_R_ptr[i] + GFPm_yp_M_R_ptr[i];
  }

  if(mlen != flen){error("lengths of male and female vectors not the same!");}

  int Total[mlen];
  for(int i=0; i<mlen; i++){
    Total[i] = Total_F[i] + Total_M[i];
  }

  /* initial genotype numbers */
  int HY[(gens+1)];
  HY[0] = 50;
  int RY[(gens+1)];
  RY[0] = 0;
  int BY[(gens+1)];
  BY[0] = 0;
  int WY[(gens+1)];
  WY[0] = 0;
  int HH[(gens+1)];
  HH[0] = 0;
  int HR[(gens+1)];
  HR[0] = 0;
  int HB[(gens+1)];
  HB[0] = 0;
  int HW[(gens+1)];
  HW[0] = 0;
  int RR[(gens+1)];
  RR[0] = 0;
  int RB[(gens+1)];
  RB[0] = 0;
  int RW[(gens+1)];
  RW[0] = 0;
  int BB[(gens+1)];
  BB[0] = 0;
  int BW[(gens+1)];
  BW[0] = 0;
  int WW[(gens+1)];
  WW[0] = 50;

  /* predicted phenotype numbers */
  int GFPp_ym_F_Pred = HH[0] + HB[0] + HW[0];
  int GFPp_yp_F_Pred = HR[0];
  int GFPm_ym_F_Pred = BB[0];
  int GFPm_yp_F_Pred = RR[0] + RB[0] + RW[0] + BW[0] + WW[0];

  int GFPp_ym_M_Pred = HY[0];
  int GFPp_yp_M_Pred = 0;
  int GFPm_ym_M_Pred = BY[0];
  int GFPm_yp_M_Pred = RY[0] + WY[0];

  /* poulation dynamic model */
  for(int i=1; i < (gens+1); i++){

  }



  /* loglikelihood of input parameters */
  double loglike = 0.;

  return ScalarReal(loglike);
};
