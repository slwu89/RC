/* ##################################################
#   MCR (Mutagenic Chain Reaction) Metropolis-Hastings MCMC
#   Sean Wu
################################################## */

#include "MCR_MetropolisHastings.h"
#include "MCR.h"
#include "rv_truncnorm.h"

/* progress bar */
void printProgress (double percentage){
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}

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


/* Random Walk Metropolis-Hastings MCMC */
SEXP C_mcmcMH(SEXP initTheta_R, SEXP proposalSD_R, SEXP numIterations_R,
              SEXP GFPp_ym_F_1_R, SEXP GFPp_yp_F_1_R, SEXP GFPm_ym_F_1_R, SEXP GFPm_yp_F_1_R,
              SEXP GFPp_ym_M_1_R, SEXP GFPp_yp_M_1_R, SEXP GFPm_ym_M_1_R, SEXP GFPm_yp_M_1_R,
              SEXP GFPp_ym_F_2_R, SEXP GFPp_yp_F_2_R, SEXP GFPm_ym_F_2_R, SEXP GFPm_yp_F_2_R,
              SEXP GFPp_ym_M_2_R, SEXP GFPp_yp_M_2_R, SEXP GFPm_ym_M_2_R, SEXP GFPm_yp_M_2_R,
              SEXP GFPp_ym_F_3_R, SEXP GFPp_yp_F_3_R, SEXP GFPm_ym_F_3_R, SEXP GFPm_yp_F_3_R,
              SEXP GFPp_ym_M_3_R, SEXP GFPp_yp_M_3_R, SEXP GFPm_ym_M_3_R, SEXP GFPm_yp_M_3_R,
              SEXP GFPp_ym_F_4_R, SEXP GFPp_yp_F_4_R, SEXP GFPm_ym_F_4_R, SEXP GFPm_yp_F_4_R,
              SEXP GFPp_ym_M_4_R, SEXP GFPp_yp_M_4_R, SEXP GFPm_ym_M_4_R, SEXP GFPm_yp_M_4_R,
              SEXP gens_R){

  /* use R's PRNG */
  GetRNGstate();

  /* calls to PROTECT */
  int protectCalls = 0;

  /* Evaluate the log posterior at initTheta: */
  double logPosteriorThetaCurrent = asReal(C_logPosterior(ScalarReal(REAL_ELT(initTheta_R,0)),
                                                   ScalarReal(REAL_ELT(initTheta_R,1)),
                                                   ScalarReal(REAL_ELT(initTheta_R,2)),
                                                   ScalarReal(REAL_ELT(initTheta_R,3)),
                                                   ScalarReal(REAL_ELT(initTheta_R,4)),
                                                   ScalarReal(REAL_ELT(initTheta_R,5)),
                                                   GFPp_ym_F_1_R, GFPp_yp_F_1_R, GFPm_ym_F_1_R, GFPm_yp_F_1_R,
                                                   GFPp_ym_M_1_R, GFPp_yp_M_1_R, GFPm_ym_M_1_R, GFPm_yp_M_1_R,
                                                   GFPp_ym_F_2_R, GFPp_yp_F_2_R, GFPm_ym_F_2_R, GFPm_yp_F_2_R,
                                                   GFPp_ym_M_2_R, GFPp_yp_M_2_R, GFPm_ym_M_2_R, GFPm_yp_M_2_R,
                                                   GFPp_ym_F_3_R, GFPp_yp_F_3_R, GFPm_ym_F_3_R, GFPm_yp_F_3_R,
                                                   GFPp_ym_M_3_R, GFPp_yp_M_3_R, GFPm_ym_M_3_R, GFPm_yp_M_3_R,
                                                   GFPp_ym_F_4_R, GFPp_yp_F_4_R, GFPm_ym_F_4_R, GFPm_yp_F_4_R,
                                                   GFPp_ym_M_4_R, GFPp_yp_M_4_R, GFPm_ym_M_4_R, GFPm_yp_M_4_R,
                                                   gens_R));

  // printf("logPosteriorThetaCurrent: %f \n",logPosteriorThetaCurrent);

  /* number of iterations to run the markov chain */
  int nIter = asInteger(numIterations_R);

  // printf("nIter %i \n",nIter);

  /* initialise variables */
  SEXP c_current = PROTECT(ScalarReal(REAL_ELT(initTheta_R,0)));
  SEXP pH_current = PROTECT(ScalarReal(REAL_ELT(initTheta_R,1)));
  SEXP prR_current = PROTECT(ScalarReal(REAL_ELT(initTheta_R,2)));
  SEXP sH_current = PROTECT(ScalarReal(REAL_ELT(initTheta_R,3)));
  SEXP sR_current = PROTECT(ScalarReal(REAL_ELT(initTheta_R,4)));
  SEXP sB_current = PROTECT(ScalarReal(REAL_ELT(initTheta_R,5)));
  protectCalls += 6;
  SEXP samples = PROTECT(allocMatrix(REALSXP,nIter,6));
  protectCalls += 1;
  int accepted = 0;

  /* proposed parameters */
  SEXP c_proposed = PROTECT(ScalarReal(0.0));
  SEXP pH_proposed = PROTECT(ScalarReal(0.0));
  SEXP prR_proposed = PROTECT(ScalarReal(0.0));
  SEXP sH_proposed = PROTECT(ScalarReal(0.0));
  SEXP sR_proposed = PROTECT(ScalarReal(0.0));
  SEXP sB_proposed = PROTECT(ScalarReal(0.0));
  protectCalls += 6;

  /* parameters of proposal distribution */
  double c_sd = REAL_ELT(proposalSD_R,0);
  double pH_sd = REAL_ELT(proposalSD_R,1);
  double prR_sd = REAL_ELT(proposalSD_R,2);
  double sH_sd = REAL_ELT(proposalSD_R,3);
  double sR_sd = REAL_ELT(proposalSD_R,4);
  double sB_sd = REAL_ELT(proposalSD_R,5);

  /* Run the MCMC algorithm for numIterations interations: */
  for(int i=0; i < nIter; i++){

    /* draw a new theta from a Gaussian proposal distribution and
     * assign this to a variable called thetaProposed
     */
    c_proposed = ScalarReal(r_truncnorm(0,1,asReal(c_current),c_sd));
    pH_proposed = ScalarReal(r_truncnorm(0,1,asReal(pH_current),pH_sd));
    prR_proposed = ScalarReal(r_truncnorm(0,1,asReal(prR_current),prR_sd));
    sH_proposed = ScalarReal(r_truncnorm(-2,1,asReal(sH_current),sH_sd));
    sR_proposed = ScalarReal(r_truncnorm(-2,1,asReal(sR_current),sR_sd));
    sB_proposed = ScalarReal(r_truncnorm(-2,1,asReal(sB_current),sB_sd));

    // printf("c_proposed: %f \n",asReal(c_proposed));

    /* Evaluate the log posterior function at the proposed theta value: */
    double logPosteriorThetaProposed = asReal(C_logPosterior(c_proposed,pH_proposed,prR_proposed,sH_proposed,sR_proposed,sB_proposed,
                                                      GFPp_ym_F_1_R, GFPp_yp_F_1_R, GFPm_ym_F_1_R, GFPm_yp_F_1_R,
                                                      GFPp_ym_M_1_R, GFPp_yp_M_1_R, GFPm_ym_M_1_R, GFPm_yp_M_1_R,
                                                      GFPp_ym_F_2_R, GFPp_yp_F_2_R, GFPm_ym_F_2_R, GFPm_yp_F_2_R,
                                                      GFPp_ym_M_2_R, GFPp_yp_M_2_R, GFPm_ym_M_2_R, GFPm_yp_M_2_R,
                                                      GFPp_ym_F_3_R, GFPp_yp_F_3_R, GFPm_ym_F_3_R, GFPm_yp_F_3_R,
                                                      GFPp_ym_M_3_R, GFPp_yp_M_3_R, GFPm_ym_M_3_R, GFPm_yp_M_3_R,
                                                      GFPp_ym_F_4_R, GFPp_yp_F_4_R, GFPm_ym_F_4_R, GFPm_yp_F_4_R,
                                                      GFPp_ym_M_4_R, GFPp_yp_M_4_R, GFPm_ym_M_4_R, GFPm_yp_M_4_R,
                                                      gens_R));

    //  printf("logPosteriorThetaProposed: %f \n",logPosteriorThetaProposed);

    /* Compute the Metropolis-Hastings (log) acceptance probability: */
    double logAcceptance = logPosteriorThetaProposed - logPosteriorThetaCurrent;

    // printf("logAcceptance: %f \n",logAcceptance);

    /* Use a random number to determine if thetaProposed will be accepted: */
    double randNum = runif(0,1);

    // printf("randNum: %f \n",randNum);
    // printf("exp(logAcceptance): %f \n",exp(logAcceptance));

    /* If accepted, update the thetaCurrent vector, etc.: */
    if(randNum < exp(logAcceptance)){
      /* update parameters */
      c_current = c_proposed;
      pH_current = pH_proposed;
      prR_current = prR_proposed;
      sH_current = sH_proposed;
      sR_current = sR_proposed;
      sB_current = sB_proposed;
      /* update target log posterior */
      logPosteriorThetaCurrent = logPosteriorThetaProposed;
      /* acceptance */
      accepted += 1;
    }

    /* Add the current theta to the vector of samples: */
    SET_REAL_ELT(samples,i+0*nIter,asReal(c_current));
    SET_REAL_ELT(samples,i+1*nIter,asReal(pH_current));
    SET_REAL_ELT(samples,i+2*nIter,asReal(prR_current));
    SET_REAL_ELT(samples,i+3*nIter,asReal(sH_current));
    SET_REAL_ELT(samples,i+4*nIter,asReal(sR_current));
    SET_REAL_ELT(samples,i+5*nIter,asReal(sB_current));

    /* progress bar */
    if(i % 100 == 0){
      printProgress((double)i/nIter);
    }
  }

  printProgress(1.0);

  /* stop using R's PRNG */
  PutRNGstate();

  /* prepare output */
  SEXP names = PROTECT(Rf_allocVector(STRSXP,2));
  protectCalls += 1;
  SET_STRING_ELT(names,0,mkChar("samples"));
  SET_STRING_ELT(names,1,mkChar("acceptance"));

  SEXP out = PROTECT(allocVector(VECSXP, 2));
  protectCalls += 1;
  namesgets(out,names);

  SET_VECTOR_ELT(out,0,samples);
  double acc = (double)accepted/nIter;
  SET_VECTOR_ELT(out,1,ScalarReal(acc));

  /* unprotect SEXPs for gc */
  UNPROTECT(protectCalls);

  /* return to R */
  return out;
};
