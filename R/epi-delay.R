##################################################
#   MGDrivE-Epi (Secret Projects)
#   fixed discrete delay
#   Sean Wu
##################################################


#' testing a fn
#'
#' this function's source is actually in utilities.h/c
#'
#' @examples
#' \dontrun{
#'
#' }
#' @useDynLib RC oneDay_oviposit_C
#' @export
oneDay_oviposit <- function(cube){
  .Call(oneDay_oviposit_C,cube)
}
