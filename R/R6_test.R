##################################################
#   Testing C code for R6 objects
#   Sean Wu
##################################################

#' Look up a symbol in an environment and print it
#'
#' @useDynLib RC find_var_C
#' @export
find_var <- function(env,var){
  invisible(.Call(find_var_C,env,var,new.env()))
}
