\name{GERGM-package}
\alias{GERGM-package}
\alias{GERGM}
\docType{package}
\title{
\packageTitle{GERGM}
}
\description{
\packageDescription{GERGM}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{GERGM}
\packageIndices{GERGM}
To use this package, first load in the network you wish to use as a (square) matrix, following the example provided below. You may then use the gergm() function to estimate a model using any combination of the following statistics: "out2star", "in2star", "ctriads", "recip", "ttriads",
"edgeweight". The gergm() function will provide all of the estimation and diagnostic functionality and the parameters of this function can be querried by typing ?gergm into the R console. Covariates are currently incorporated into the model but have not been documented or thoroughly tested as of 6/23/15. This will be corrected shortly. 
}
\author{
\packageAuthor{GERGM}
Maintainer: \packageMaintainer{GERGM}
}
\references{
Relevant papers detailing the model can be found at the links below. Please cite the second paper when using this package.

* Bruce A. Desmarais, and Skyler J. Cranmer,  (2012). "Statistical inference for valued-edge networks: the generalized exponential random graph model". PloS One. http://dx.plos.org/10.1371/journal.pone.0030136)

* James D. Wilson, Matthew J. Denny, Shankar Bhamidi, Skyler Cranmer, and Bruce Desmarais (2015). "Stochastic Weighted Graphs: Flexible Model Specification and Simulation". http://arxiv.org/abs/1505.04015
}
\examples{
net <- matrix(runif(100),10,10)
diag(net) <- 0
formula <- "net ~ recip + edgeweight"  
      
test <- gergm(formula,
              network_is_directed = TRUE,
              use_MPLE_only = FALSE,
              data_transformation = NULL,
              estimation_method = "Gibbs",
              maximum_number_of_lambda_updates = 1,
              maximum_number_of_theta_updates = 5,
              number_of_networks_to_simulate = 5000,
              thin = 1/10,
              proposal_variance = 0.5,
              exponential_weights = NULL,
              downweight_statistics_together = TRUE,
              MCMC_burnin = 1000,
              seed = 123,
              convergence_tolerance = 0.01,
              MPLE_gain_factor = 0,
              force_x_theta_update = 2)
}