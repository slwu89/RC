##################################################
#   interacting with environments tidbits
#   Sean Wu
##################################################

#' apply a function f to all elements in a hashed envir for side-effects (no return value)
#'
#' @examples
#' \dontrun{
#' e <- new.env()
#' for(a in LETTERS[1:10]){
#'   assign(x = a,value = rpois(n = 5,lambda = 10),envir = e)
#'  }
#'  fn <- function(x,msg){
#'  cat(msg,"\n")
#'  cat(x,"\n")
#'  }
#' eapply_fast(env = e,fn = fn,msg = "hi")
#' }
#' @useDynLib RC eapply_fast_C
#' @export
eapply_fast <- function(env,fn,...){
  invisible(.Call(eapply_fast_C,env,fn,new.env()))
}
