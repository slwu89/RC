#' MCR model likelihood (interface to C)
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


#' MCR model likelihood for all experiments (interface to C)
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
