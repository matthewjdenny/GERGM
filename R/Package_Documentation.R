#' GERGM: Generalized Exponential Random Graph Model
#'
#' @section GERGM functions:
#' To use this package, first load in the network you wish to use as a (square)
#' matrix, following the example provided below. You may then use the gergm()
#' function to estimate a model using any combination of the following statistics:
#' "out2stars", "in2stars", "ctriads", "mutual", "ttriads", "edges",
#' "absdiff(covariate)", "edgecov(covariate)", "sender(covariate)",
#' "reciever(covariate)", "nodematch(covariate)", "nodemix(covariate)",
#' "netcov(network_covariate)". The gergm() function provides all of the basic
#' estimation and diagnostic functionality and the parameters of this function
#' can be queried by typing ?gergm into the R console. If you wish to access
#' additional fit and degeneracy diagnostic functionality, the GOF(),
#' Estimate_Plot(), Trace_Plot() and hysteresis() functions can be accessed
#' directly by the user.  You may also plot the initial network using
#' plot_network() and simulate networks for given structural parameters using
#' the simulate_networks() function. Experimental support for specifying
#' multiple GERGMs in parallel (allowing for different equations, dependent
#' networks and covariates) is available in the parallel_gergm() function.
#' An experimental feature which seeks to automatically optimize model
#' hyperparameters for best fit and to attempt to deal with degeneracy issues
#' can be turned on be specifying hyperparameter_optimization = TRUE.
#'
#' @docType package
#' @name GERGM
NULL
#> NULL

#' @import methods
NULL

#' @importFrom grDevices dev.off gray pdf rgb colorRampPalette
NULL

#' @importFrom graphics boxplot legend lines par plot text axis plot.new layout abline points
NULL

#' @import plyr
NULL

#' @importFrom  stats as.formula dgamma dt lm optim pgamma pnorm pt qgamma qnorm qt rnorm runif sd t.test cor dbeta pbeta qbeta density
NULL

#' @importFrom utils combn write.table
NULL

#' @useDynLib GERGM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL
#> NULL
