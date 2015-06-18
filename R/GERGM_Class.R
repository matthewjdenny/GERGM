# an S4 class for gergm objects
setClass(Class = "gergm", 
         representation = representation(
           network = "matrix",
           bounded.network = "matrix",
           formula = "formula",
           stats = "matrix",
           theta.coef = "data.frame", 
           lambda.coef = "data.frame",
           weights = "numeric"
         ), 
         validity = function(object) {
           if (!"matrix" %in% class(object@network) & is.null(object@network) 
               == FALSE) {
             stop("'network' must be a 'matrix' object or 'NULL'.")
           }
           if (!"matrix" %in% class(object@bounded.network)) {
             stop("'bounded.network' must be a 'matrix' object.")
           }
           if (!is.data.frame(object@theta.coef)) {
             stop("'theta.coef' must be a 'data.frame' object.")
           }
           if (!"formula" %in% class(object@formula)) {
             stop("'formula' is not a 'formula' object.")
           }
           if (!is.data.frame(object@lambda.coef) & is.null(object@lambda.coef) 
               == FALSE) {
             stop("'lambda.coef' must be a 'data.frame' object or 'NULL'.")
           }
           if (!"matrix" %in% class(object@stats)) {
             stop("'stats' must be a 'matrix' object.")
           }
           if (!is.numeric(object@weights) & is.null(object@weights) 
               == FALSE) {
             stop("'weights' must be a 'data.frame' object or 'NULL'.")
           }
           return(TRUE)
         }
)

# define coef for pretty output of gergm object
setMethod(f = "coef", signature = "gergm", definition = function(object, ...) {
  return(list(Theta  = object@theta.coef, Lambda = object@lambda.coef))
}
)

# define 'show' to get a pretty output of the gergm
setMethod(f = "show", signature = "gergm", definition = function(object){
  message("Theta:")
  print(object@theta.coef)
  message("Lambda:")
  print(object@lambda.coef)
  message("Weights:")
  print(object@weights)
  message("Network Statistics:")
  print(object@stats)
}
)