# an S4 class for gergm objects
setClass(Class = "gergm",
         representation = representation(
           network = "matrix",
           bounded.network = "matrix",
           formula = "formula",
           stats = "matrix",
           theta.coef = "data.frame",
           lambda.coef = "data.frame",
           weights = "numeric",
           num_nodes = "numeric",
           MCMC_output = "list",
           observed_network = "matrix",
           observed_bounded_network = "matrix",
           data_transformation = "array"
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
           if (!is.numeric(object@num_nodes) & is.null(object@num_nodes)
               == FALSE) {
             stop("'num_nodes' must be a 'data.frame' object or 'NULL'.")
           }
           if (!"list" %in% class(object@MCMC_output)) {
             stop("'MCMC_output' must be a 'list' object.")
           }
           if (!"matrix" %in% class(object@observed_network)) {
             stop("'observed_network' must be a 'matrix' object.")
           }
           if (!"matrix" %in% class(object@observed_bounded_network)) {
             stop("'observed_bounded_network' must be a 'matrix' object.")
           }
           if (!"array" %in% class(object@data_transformation)) {
             stop("'data_transformation' must be a 'array' object.")
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
