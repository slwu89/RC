##################################################
#   MCR (Mutagenic Chain Reaction) parameter inference
#   Sean Wu
##################################################

#' MCR model log likelihood (interface to C)
#'
#' @useDynLib RC C_logLike_MCRMod
#' @export
logLike_MCRMod <- function(c, pH, prR, sH, sR, sB,
                           GFPp_ym_F, GFPp_yp_F, GFPm_ym_F, GFPm_yp_F,
                           GFPp_ym_M, GFPp_yp_M, GFPm_ym_M, GFPm_yp_M,
                           gens){
  .Call(C_logLike_MCRMod,as.double(c),as.double(pH),as.double(prR),as.double(sH),as.double(sR),as.double(sB),
         as.double(GFPp_ym_F),as.double(GFPp_yp_F),as.double(GFPm_ym_F),as.double(GFPm_yp_F),
         as.double(GFPp_ym_M),as.double(GFPp_yp_M),as.double(GFPm_ym_M),as.double(GFPm_yp_M),
         as.double(gens))
}


#' MCR model log likelihood for all experiments (interface to C)
#'
#' @useDynLib RC C_logLike_AllExpts
#' @export
logLike_AllExpts <- function(c, pH, prR, sH, sR, sB,
                           GFPp_ym_1_F, GFPp_yp_1_F, GFPm_ym_1_F, GFPm_yp_1_F,
                           GFPp_ym_1_M, GFPp_yp_1_M, GFPm_ym_1_M, GFPm_yp_1_M,
                           GFPp_ym_2_F, GFPp_yp_2_F, GFPm_ym_2_F, GFPm_yp_2_F,
                           GFPp_ym_2_M, GFPp_yp_2_M, GFPm_ym_2_M, GFPm_yp_2_M,
                           GFPp_ym_3_F, GFPp_yp_3_F, GFPm_ym_3_F, GFPm_yp_3_F,
                           GFPp_ym_3_M, GFPp_yp_3_M, GFPm_ym_3_M, GFPm_yp_3_M,
                           GFPp_ym_4_F, GFPp_yp_4_F, GFPm_ym_4_F, GFPm_yp_4_F,
                           GFPp_ym_4_M, GFPp_yp_4_M, GFPm_ym_4_M, GFPm_yp_4_M,
                           gens){
  .Call(C_logLike_AllExpts,as.double(c),as.double(pH),as.double(prR),as.double(sH),as.double(sR),as.double(sB),
        as.double(GFPp_ym_1_F),as.double(GFPp_yp_1_F),as.double(GFPm_ym_1_F),as.double(GFPm_yp_1_F),
        as.double(GFPp_ym_1_M),as.double(GFPp_yp_1_M),as.double(GFPm_ym_1_M),as.double(GFPm_yp_1_M),
        as.double(GFPp_ym_2_F),as.double(GFPp_yp_2_F),as.double(GFPm_ym_2_F),as.double(GFPm_yp_2_F),
        as.double(GFPp_ym_2_M),as.double(GFPp_yp_2_M),as.double(GFPm_ym_2_M),as.double(GFPm_yp_2_M),
        as.double(GFPp_ym_3_F),as.double(GFPp_yp_3_F),as.double(GFPm_ym_3_F),as.double(GFPm_yp_3_F),
        as.double(GFPp_ym_3_M),as.double(GFPp_yp_3_M),as.double(GFPm_ym_3_M),as.double(GFPm_yp_3_M),
        as.double(GFPp_ym_4_F),as.double(GFPp_yp_4_F),as.double(GFPm_ym_4_F),as.double(GFPm_yp_4_F),
        as.double(GFPp_ym_4_M),as.double(GFPp_yp_4_M),as.double(GFPm_ym_4_M),as.double(GFPm_yp_4_M),
        as.double(gens))
}

#' MCR model log prior distribution (interface to C)
#'
#' @useDynLib RC C_logPrior
#' @export
logPrior <- function(c, pH, prR,sH, sR, sB){
  .Call(C_logPrior,as.double(c),as.double(pH),as.double(prR),as.double(sH),as.double(sR),as.double(sB))
}

#' MCR model log target posterior distribution (interface to C)
#'
#' @useDynLib RC C_logPosterior
#' @export
logPosterior <- function(c, pH, prR, sH, sR, sB,
                         GFPp_ym_F_1, GFPp_yp_F_1, GFPm_ym_F_1, GFPm_yp_F_1,
                         GFPp_ym_M_1, GFPp_yp_M_1, GFPm_ym_M_1, GFPm_yp_M_1,
                         GFPp_ym_F_2, GFPp_yp_F_2, GFPm_ym_F_2, GFPm_yp_F_2,
                         GFPp_ym_M_2, GFPp_yp_M_2, GFPm_ym_M_2, GFPm_yp_M_2,
                         GFPp_ym_F_3, GFPp_yp_F_3, GFPm_ym_F_3, GFPm_yp_F_3,
                         GFPp_ym_M_3, GFPp_yp_M_3, GFPm_ym_M_3, GFPm_yp_M_3,
                         GFPp_ym_F_4, GFPp_yp_F_4, GFPm_ym_F_4, GFPm_yp_F_4,
                         GFPp_ym_M_4, GFPp_yp_M_4, GFPm_ym_M_4, GFPm_yp_M_4,
                         gens){
  .Call(C_logPosterior,as.double(c),as.double(pH),as.double(prR),as.double(sH),as.double(sR),as.double(sB),
        as.double(GFPp_ym_1_F),as.double(GFPp_yp_1_F),as.double(GFPm_ym_1_F),as.double(GFPm_yp_1_F),
        as.double(GFPp_ym_1_M),as.double(GFPp_yp_1_M),as.double(GFPm_ym_1_M),as.double(GFPm_yp_1_M),
        as.double(GFPp_ym_2_F),as.double(GFPp_yp_2_F),as.double(GFPm_ym_2_F),as.double(GFPm_yp_2_F),
        as.double(GFPp_ym_2_M),as.double(GFPp_yp_2_M),as.double(GFPm_ym_2_M),as.double(GFPm_yp_2_M),
        as.double(GFPp_ym_3_F),as.double(GFPp_yp_3_F),as.double(GFPm_ym_3_F),as.double(GFPm_yp_3_F),
        as.double(GFPp_ym_3_M),as.double(GFPp_yp_3_M),as.double(GFPm_ym_3_M),as.double(GFPm_yp_3_M),
        as.double(GFPp_ym_4_F),as.double(GFPp_yp_4_F),as.double(GFPm_ym_4_F),as.double(GFPm_yp_4_F),
        as.double(GFPp_ym_4_M),as.double(GFPp_yp_4_M),as.double(GFPm_ym_4_M),as.double(GFPm_yp_4_M),
        as.double(gens))
}
