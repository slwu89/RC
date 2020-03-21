##################################################
#   testing how to copy subsets of matrices
#   Sean Wu
##################################################

#' Subset Matrix in C by rows
#'
#'
#' @examples
#' \dontrun{
#'    mat <- matrix(1:15,nrow=5,ncol=3)
#'    submat(mat,2)
#' }
#' @useDynLib RC submat_C
#' @export
submat <- function(mat,nrow){
  storage.mode(mat) <- "double"
  .Call(submat_C,as.matrix(mat),as.integer(nrow),as.integer(ncol(mat)))
}

#' Subset Matrix in C by cols
#'
#'
#' @examples
#' \dontrun{
#'    mat <- matrix(1:15,nrow=3,ncol=5)
#'    submat_cols(mat,2)
#' }
#' @useDynLib RC submat_cols_C
#' @export
submat_cols <- function(mat,ncol){
  storage.mode(mat) <- "double"
  .Call(submat_cols_C,as.matrix(mat),as.integer(nrow(mat)),as.integer(ncol))
}
