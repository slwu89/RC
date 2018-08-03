#' MCR model likelihood (interface to C)
#'
#' @useDynLib RC C_logLike_MCRMod
#' @export
logLike_MCRMod <- function(c, pH, prR, sH, sR, sB,
                           GFPp_ym_F, GFPp_yp_F, GFPm_ym_F, GFPm_yp_F,
                           GFPp_ym_M, GFPp_yp_M, GFPm_ym_M, GFPm_yp_M,
                           gens){
  .Call(C_logLike_MCRMod,as.double(c),as.double(pH),as.double(prR),as.double(sH),as.double(sR),as.double(sB),
         as.integer(GFPp_ym_F),as.integer(GFPp_yp_F),as.integer(GFPm_ym_F),as.integer(GFPm_yp_F),
         as.integer(GFPp_ym_M),as.integer(GFPp_yp_M),as.integer(GFPm_ym_M),as.integer(GFPm_yp_M),
         as.integer(gens))
}
