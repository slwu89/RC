##################################################
#   Functional Programming tidbits
#   Sean Wu
##################################################

#' Simple version of "Reduce"
#'
#' simple version of Reduce (see http://adv-r.had.co.nz/Functionals.html#functionals-fp)
#'
#' @examples
#' \dontrun{
#'
#' }
#' @useDynLib RC Reduce_Simple_C
#' @export
Reduce_Simple <- function(f,x){
  f <- match.fun(f)
  .Call(Reduce_Simple_C,f,x,new.env())
}
