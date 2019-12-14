

#' MCSM: Algorithm 1 "Pool Adjacent Violators"
#'
#' blah
#'
#' @examples
#' \dontrun{
#' f <- c(23,27,25,28)
#' w <- rep(1/4,4)
#'
#' }
#' @useDynLib RC mcsm_pava_c
#' @export
mcsm_pava <- function(f,w){
  .Call(mcsm_pava_c,f,w)
}
