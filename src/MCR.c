#include "MCR.h"

/* MCR model likelihood */
SEXP C_logLike_MCRMod(SEXP c_R, SEXP pH_R, SEXP prR_R, SEXP sH_R, SEXP sR_R, SEXP sB_R,
                      SEXP GFPp_ym_F_R, SEXP GFPp_yp_F_R, SEXP GFPm_ym_F_R, SEXP GFPm_yp_F_R,
                      SEXP GFPp_ym_M_R, SEXP GFPp_yp_M_R, SEXP GFPm_ym_M_R, SEXP GFPm_yp_M_R)
{

  /* model parameters as C types */
  double c = asReal(c_R);
  double pH = asReal(pH_R);
  double prR = asReal(prR_R);
  double sH = asReal(sH_R);
  double sR = asReal(sR_R);
  double sB = asReal(sB_R);

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
  int HY = 50;
  int RY = 0;
  int BY = 0;
  int WY = 0;
  int HH = 0;
  int HR = 0;
  int HB = 0;
  int HW = 0;
  int RR = 0;
  int RB = 0;
  int RW = 0;
  int BB = 0;
  int BW = 0;
  int WW = 50;

  /* predicted phenotype numbers */
  int GFPp_ym_F_Pred = HH + HB + HW;
  int GFPp_yp_F_Pred = HR;
  int GFPm_ym_F_Pred = BB;
  int GFPm_yp_F_Pred = RR + RB + RW + BW + WW;

  int GFPp_ym_M_Pred = HY;
  int GFPp_yp_M_Pred = 0;
  int GFPm_ym_M_Pred = BY;
  int GFPm_yp_M_Pred = RY + WY;


  /* loglikelihood of input parameters */
  double loglike = 0.;

  return ScalarReal(loglike);
};
