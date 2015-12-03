#' GERGM: A package for estimating and diagnosing Generalized Exponential Random
#' Graph Models
#'
#' @section GERGM functions:
#' To use this package, first load in the network you wish to use as a (square)
#' matrix, following the example provided below. You may then use the gergm()
#' function to estimate a model using any combination of the following statistics:
#' "out2star", "in2star", "ctriads", "recip", "ttriads", "edges",
#' "absdiff(covariate)", "edgecov(covariate)", "sender(covariate)",
#' "reciever(covariate)", "nodefactor(covariate)", "netcov(network_covariate)".
#' The gergm() function will provide all of the estimation and diagnostic
#' functionality and the parameters of this function can be querried by typing
#' ?gergm into the R console. You may also plot the initial network using
#' plot_network() and simulate networks using the simulate_networks() function.
#'
#' @docType package
#' @name GERGM
NULL
#> NULL

#' @import methods
#' @import RcppArmadillo
NULL

#' @importFrom grDevices dev.off gray pdf rgb colorRampPalette
NULL

#' @importFrom graphics boxplot legend lines par plot text
NULL

#' @importFrom  stats as.formula dgamma dt lm optim pgamma pnorm pt qgamma qnorm qt rnorm runif sd t.test
NULL

#' @importFrom utils combn write.table
NULL

#' @useDynLib GERGM
#' @importFrom Rcpp sourceCpp
NULL
#> NULL
