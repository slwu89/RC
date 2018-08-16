##################################################
#   Testing C code for R6 objects
#   Sean Wu
##################################################

R6testClassGenerator <- R6:::R6Class("R6testClassGenerator",
                  public = list(
                    initialize = function(param){
                      cat("creating R6testClassGenerator at: ",pryr:::address(self),"\n")
                      private$param = param
                    },
                    finalize = function(){
                      cat("destroying R6testClassGenerator at: ",pryr:::address(self),"\n")
                    },
                    get_param = function(){return(private$param)}
                  ),
                  private = list(
                    param = numeric(1)
                  )
)

R6testClass <- R6testClassGenerator$new(5)

#' Look up a symbol in an environment and print it
#'
#' @examples
#' \dontrun{
#' library(R6)
#' Person <- R6Class("Person",
#'                   public = list(
#'                     name = NULL,
#'                     hair = NULL,
#'                     initialize = function(name = NA, hair = NA) {
#'                       self$name <- name
#'                       self$hair <- hair
#'                       self$greet()
#'                     },
#'                     set_hair = function(val) {
#'                       self$hair <- val
#'                     },
#'                     greet = function() {
#'                       cat(paste0("Hello, my name is ", self$name, ".\n"))
#'                     },
#'                     greet2 = function(){
#'                       RC::find_var(private,"z")
#'                     }
#'                   ),
#'                   private = list(z=5)
#' )
#'
#' x <- Person$new("bob","black")
#' x$greet2()
#' }
#' @useDynLib RC find_var_C
#' @export
find_var <- function(env,var){
  invisible(.Call(find_var_C,env,var,new.env()))
}
