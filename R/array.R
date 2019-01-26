##################################################
#   interacting with arrays tidbits
#   Sean Wu
##################################################

#' subset a 3d array
#'
#' @examples
#' \dontrun{
#'  x <- array(data = 1:(2*3*4),dim = c(2,3,4))
#' x[1,3,4]
#' sub_array3(1,3,4,x)
#' }
#'
#' @useDynLib RC sub_array3_C
#' @export
sub_array3 <- function(i,j,k,array){
  .Call(sub_array3_C,as.integer(i),as.integer(j),as.integer(k),array)
}

#' subset a 4d array
#'
#' @examples
#' \dontrun{
#'  x <- array(data = 1:(2*3*4*5),dim = c(2,3,4,5))
#' x[1,3,4,5]
#' sub_array4(1,3,4,5,x)
#' }
#'
#' @useDynLib RC sub_array4_C
#' @export
sub_array4 <- function(i,j,k,l,array){
  .Call(sub_array4_C,as.integer(i),as.integer(j),as.integer(k),as.integer(l),array)
}
