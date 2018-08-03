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
  double HY[(gens+1)];
  HY[0] = 50;
  double RY[(gens+1)];
  RY[0] = 0;
  double BY[(gens+1)];
  BY[0] = 0;
  double WY[(gens+1)];
  WY[0] = 0;
  double HH[(gens+1)];
  HH[0] = 0;
  double HR[(gens+1)];
  HR[0] = 0;
  double HB[(gens+1)];
  HB[0] = 0;
  double HW[(gens+1)];
  HW[0] = 0;
  double RR[(gens+1)];
  RR[0] = 0;
  double RB[(gens+1)];
  RB[0] = 0;
  double RW[(gens+1)];
  RW[0] = 0;
  double BB[(gens+1)];
  BB[0] = 0;
  double BW[(gens+1)];
  BW[0] = 0;
  double WW[(gens+1)];
  WW[0] = 50;

  /* predicted phenotype numbers */
  double GFPp_ym_F_Pred[(gens+1)];
  GFPp_ym_F_Pred[0] = HH[0] + HB[0] + HW[0];
  double GFPp_yp_F_Pred[(gens+1)];
  GFPp_yp_F_Pred[0] = HR[0];
  double GFPm_ym_F_Pred[(gens+1)];
  GFPm_ym_F_Pred[0] = BB[0];
  double GFPm_yp_F_Pred[(gens+1)];
  GFPm_yp_F_Pred[0] = RR[0] + RB[0] + RW[0] + BW[0] + WW[0];

  double GFPp_ym_M_Pred[(gens+1)];
  GFPp_ym_M_Pred[0] = HY[0];
  double GFPp_yp_M_Pred[(gens+1)];
  GFPp_yp_M_Pred[0] = 0;
  double GFPm_ym_M_Pred[(gens+1)];
  GFPm_ym_M_Pred[0] = BY[0];
  double GFPm_yp_M_Pred[(gens+1)];
  GFPm_yp_M_Pred[0] = RY[0] + WY[0];

  /* poulation dynamic model */
  for(int i=1; i < (gens+1); i++){

    /* daily updating populations */
    int HY_Temp = 0.5*HY[i-1]*HH[i-1] + 0.5*RY[i-1]*HH[i-1] + 0.5*BY[i-1]*HH[i-1] + 0.5*WY[i-1]*HH[i-1] +
      0.25*HY[i-1]*HR[i-1] + 0.25*RY[i-1]*HR[i-1] + 0.25*BY[i-1]*HR[i-1] + 0.25*WY[i-1]*HR[i-1] +
      0.25*HY[i-1]*HB[i-1] + 0.25*RY[i-1]*HB[i-1] + 0.25*BY[i-1]*HB[i-1] + 0.25*WY[i-1]*HB[i-1] +
      0.25*(1 + c*pH)*HY[i-1]*HW[i-1] + 0.25*(1 + c*pH)*RY[i-1]*HW[i-1] + 0.25*(1 + c*pH)*BY[i-1]*HW[i-1] + 0.25*(1 + c*pH)*WY[i-1]*HW[i-1];

    int RY_Temp = 0.5*HY[i-1]*RR[i-1] + 0.5*RY[i-1]*RR[i-1] + 0.5*BY[i-1]*RR[i-1] + 0.5*WY[i-1]*RR[i-1] +
      0.25*HY[i-1]*HR[i-1] + 0.25*RY[i-1]*HR[i-1] + 0.25*BY[i-1]*HR[i-1] + 0.25*WY[i-1]*HR[i-1] +
      0.25*HY[i-1]*RB[i-1] + 0.25*RY[i-1]*RB[i-1] + 0.25*BY[i-1]*RB[i-1] + 0.25*WY[i-1]*RB[i-1] +
      0.25*HY[i-1]*RW[i-1] + 0.25*RY[i-1]*RW[i-1] + 0.25*BY[i-1]*RW[i-1] + 0.25*WY[i-1]*RW[i-1] +
      0.25*c*pR*HY[i-1]*HW[i-1] + 0.25*c*pR*RY[i-1]*HW[i-1] + 0.25*c*pR*BY[i-1]*HW[i-1] + 0.25*c*pR*WY[i-1]*HW[i-1];

    int BY_Temp = 0.5*HY[i-1]*BB[i-1] + 0.5*RY[i-1]*BB[i-1] + 0.5*BY[i-1]*BB[i-1] + 0.5*WY[i-1]*BB[i-1] +
      0.25*HY[i-1]*HB[i-1] + 0.25*RY[i-1]*HB[i-1] + 0.25*BY[i-1]*HB[i-1] + 0.25*WY[i-1]*HB[i-1] +
      0.25*HY[i-1]*RB[i-1] + 0.25*RY[i-1]*RB[i-1] + 0.25*BY[i-1]*RB[i-1] + 0.25*WY[i-1]*RB[i-1] +
      0.25*HY[i-1]*BW[i-1] + 0.25*RY[i-1]*BW[i-1] + 0.25*BY[i-1]*BW[i-1] + 0.25*WY[i-1]*BW[i-1] +
      0.25*c*pB*HY[i-1]*HW[i-1] + 0.25*c*pB*RY[i-1]*HW[i-1] + 0.25*c*pB*BY[i-1]*HW[i-1] + 0.25*c*pB*WY[i-1]*HW[i-1];

    int WY_Temp = 0.5*HY[i-1]*WW[i-1] + 0.5*RY[i-1]*WW[i-1] + 0.5*BY[i-1]*WW[i-1] + 0.5*WY[i-1]*WW[i-1] +
      0.25*HY[i-1]*RW[i-1] + 0.25*RY[i-1]*RW[i-1] + 0.25*BY[i-1]*RW[i-1] + 0.25*WY[i-1]*RW[i-1] +
      0.25*HY[i-1]*BW[i-1] + 0.25*RY[i-1]*BW[i-1] + 0.25*BY[i-1]*BW[i-1] + 0.25*WY[i-1]*BW[i-1] +
      0.25*(1 - c)*HY[i-1]*HW[i-1] + 0.25*(1 - c)*RY[i-1]*HW[i-1] + 0.25*(1 - c)*BY[i-1]*HW[i-1] + 0.25*(1 - c)*WY[i-1]*HW[i-1];

    int HH_Temp = 0.5*HY[i-1]*HH[i-1] + 0.25*HY[i-1]*HR[i-1] + 0.25*HY[i-1]*HB[i-1] + 0.25*(1 + c*pH)*HY[i-1]*HW[i-1];

    int HR_Temp = 0.5*RY[i-1]*HH[i-1] + 0.25*HY[i-1]*HR[i-1] + 0.25*RY[i-1]*HR[i-1] + 0.25*RY[i-1]*HB[i-1] + 0.25*c*pR*HY[i-1]*HW[i-1] +
      0.25*(1 + c*pH)*RY[i-1]*HW[i-1] + 0.5*HY[i-1]*RR[i-1] + 0.25*HY[i-1]*RB[i-1] + 0.25*HY[i-1]*RW[i-1];

    int HB_Temp = 0.5*BY[i-1]*HH[i-1] + 0.25*HY[i-1]*HB[i-1] + 0.25*BY[i-1]*HB[i-1] + 0.25*BY[i-1]*HR[i-1] + 0.25*c*pB*HY[i-1]*HW[i-1] +
      0.25*(1 + c*pH)*BY[i-1]*HW[i-1] + 0.5*HY[i-1]*BB[i-1] + 0.25*HY[i-1]*RB[i-1] + 0.25*HY[i-1]*BW[i-1];

    int HW_Temp = 0.5*WY[i-1]*HH[i-1] + 0.25*(1 - c)*HY[i-1]*HW[i-1] + 0.25*(1 + c*pH)*WY[i-1]*HW[i-1] + 0.25*WY[i-1]*HR[i-1] +
      0.5*HY[i-1]*WW[i-1] + 0.25*HY[i-1]*RW[i-1] + 0.25*HY[i-1]*BW[i-1] + 0.25*WY[i-1]*HB[i-1];

    int RR_Temp = 0.5*RY[i-1]*RR[i-1] + 0.25*RY[i-1]*HR[i-1] + 0.25*RY[i-1]*RB[i-1] + 0.25*RY[i-1]*RW[i-1] + 0.25*c*pR*RY[i-1]*HW[i-1];

    int RB_Temp = 0.5*BY[i-1]*RR[i-1] + 0.25*RY[i-1]*HB[i-1] + 0.25*BY[i-1]*HR[i-1] + 0.25*BY[i-1]*RW[i-1] + 0.25*c*pB*RY[i-1]*HW[i-1] +
      0.5*RY[i-1]*BB[i-1] + 0.25*RY[i-1]*RB[i-1] + 0.25*RY[i-1]*BW[i-1] + 0.25*BY[i-1]*RB[i-1] + 0.25*c*pR*BY[i-1]*HW[i-1];

    int RW_Temp = 0.5*WY[i-1]*RR[i-1] + 0.25*(1 - c)*RY[i-1]*HW[i-1] + 0.25*c*pR*WY[i-1]*HW[i-1] + 0.25*WY[i-1]*HR[i-1] +
      0.5*RY[i-1]*WW[i-1] + 0.25*RY[i-1]*RW[i-1] + 0.25*RY[i-1]*BW[i-1] + 0.25*WY[i-1]*RB[i-1] + 0.25*WY[i-1]*RW[i-1];

    int BB_Temp = 0.5*BY[i-1]*BB[i-1] + 0.25*BY[i-1]*HB[i-1] + 0.25*BY[i-1]*RB[i-1] + 0.25*BY[i-1]*BW[i-1] + 0.25*c*pB*BY[i-1]*HW[i-1];

    int BW_Temp = 0.5*WY[i-1]*BB[i-1] + 0.25*(1 - c)*BY[i-1]*HW[i-1] + 0.25*c*pB*WY[i-1]*HW[i-1] + 0.25*WY[i-1]*HB[i-1] +
      0.5*BY[i-1]*WW[i-1] + 0.25*BY[i-1]*BW[i-1] + 0.25*BY[i-1]*RW[i-1] + 0.25*WY[i-1]*RB[i-1] + 0.25*WY[i-1]*BW[i-1];

    int WW_Temp = 0.5*WY[i-1]*WW[i-1] + 0.25*WY[i-1]*RW[i-1] + 0.25*WY[i-1]*BW[i-1] + 0.25*(1 - c)*WY[i-1]*HW[i-1];

    int W = (1-sH)*HY_Temp + (1-sR)*RY_Temp + (1-sB)*BY_Temp + WY_Temp +
      (1-sH)*(1-sH)*HH_Temp + (1-sH)*(1-sR)*HR_Temp + (1-sH)*(1-sB)*HB_Temp +
      (1-sH)*HW_Temp + (1-sR)*(1-sR)*RR_Temp + (1-sR)*(1-sB)*RB_Temp +
      (1-sR)*RW_Temp + (1-sB)*(1-sB)*BB_Temp + (1-sB)*BW_Temp + WW_Temp;

    /* push back into population arrays */
    HY[i] = (1-sH)*HY_Temp/W;
    RY[i] = (1-sR)*RY_Temp/W;
    BY[i] = (1-sB)*BY_Temp/W;
    WY[i] = WY_Temp/W;
    HH[i] = (1-sH)*(1-sH)*HH_Temp/W;
    HR[i] = (1-sH)*(1-sR)*HR_Temp/W;
    HB[i] = (1-sH)*(1-sB)*HB_Temp/W;
    HW[i] = (1-sH)*HW_Temp/W;
    RR[i] = (1-sR)*(1-sR)*RR_Temp/W;
    RB[i] = (1-sR)*(1-sB)*RB_Temp/W;
    RW[i] = (1-sR)*RW_Temp/W;
    BB[i] = (1-sB)*(1-sB)*BB_Temp/W;
    BW[i] = (1-sB)*BW_Temp/W;
    WW[i] = WW_Temp/W;

    GFPp_ym_F_Pred[i] = HH[i] + HB[i] + HW[i];
    GFPp_yp_F_Pred[i] = HR[i];
    GFPm_ym_F_Pred[i] = BB[i];
    GFPm_yp_F_Pred[i] = RR[i] + RB[i] + RW[i] + BW[i] + WW[i];

    GFPp_ym_M_Pred[i] = HY[i];
    GFPp_yp_M_Pred[i] = 0;
    GFPm_ym_M_Pred[i] = BY[i];
    GFPm_yp_M_Pred[i] = RY[i] + WY[i];
  }

  /* multinomial likelihood calculation */
  double loglike = 0.0;

  /* iterate over length of input data */
  for(int i=0 ; i<mlen; i++){
    loglike += GFPp_ym_F_R_ptr[i]*log(fmax(GFPp_ym_F_Pred[i+1],1e-10)) +
      GFPp_yp_F_R_ptr[i]*log(fmax(GFPp_yp_F_Pred[i+1],1e-10)) +
      GFPm_ym_F_R_ptr[i]*log(fmax(GFPm_ym_F_Pred[i+1],1e-10)) +
      GFPm_yp_F_R_ptr[i]*log(fmax(GFPm_yp_F_Pred[i+1],1e-10)) +
      GFPp_ym_M_R_ptr[i]*log(fmax(GFPp_ym_M_Pred[i+1],1e-10)) +
      GFPp_yp_M_R_ptr[i]*log(fmax(GFPp_yp_M_Pred[i+1],1e-10)) +
      GFPm_ym_M_R_ptr[i]*log(fmax(GFPm_ym_M_Pred[i+1],1e-10)) +
      GFPm_yp_M_R_ptr[i]*log(fmax(GFPm_yp_M_Pred[i+1],1e-10));
  }

  return ScalarReal(loglike);
};
