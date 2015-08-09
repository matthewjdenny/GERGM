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
To use this package, first load in the network you wish to use as a (square) matrix, following the example provided below. You may then use the gergm() function to estimate a model using any combination of the following statistics: "out2star(alpha = 1)", "in2star(alpha = 1)", "ctriads(alpha = 1)", "recip(alpha = 1)", "ttriads(alpha = 1)", "edges(alpha = 1)", "absdiff(covariate = "MyCov")", "edgecov(covariate = "MyCov")", "sender(covariate = "MyCov")", "reciever(covariate = "MyCov")", "nodefactor(covariate, base = "MyBase")", "netcov(network_covariate)". To use exponential downweighting for any of the network level terms, simply specify a value for alpha less than 1. The gergm() function will provide all of the estimation and diagnostic functionality and the parameters of this function can be querried by typing ?gergm into the R console. You may also generate diagnostic plots using a GERGM Object returned by the gergm() function by using any of the following functions: Estimate_Plot(), GOF(), Trace_Plot(). See examples below.
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
\dontrun{
########################### 1. No Covariates #############################
# Preparing an unbounded network without covariates for gergm estimation #
net <- matrix(rnorm(100,0,20),10,10)
colnames(net) <- rownames(net) <- letters[1:10]
formula <- net ~ recip + edges

test <- gergm(formula,
	          normalization_type = "division",
	          network_is_directed = TRUE,
	          use_MPLE_only = FALSE,
	          estimation_method = "Metropolis",
	          maximum_number_of_lambda_updates = 1,
	          maximum_number_of_theta_updates = 5,
	          number_of_networks_to_simulate = 40000,
	          thin = 1/10,
	          proposal_variance = 0.5,
	          downweight_statistics_together = TRUE,
	          MCMC_burnin = 10000,
	          seed = 456,
	          convergence_tolerance = 0.01,
	          MPLE_gain_factor = 0,
	          force_x_theta_update = 4)

########################### 2. Covariates #############################
# Preparing an unbounded network with covariates for gergm estimation #
net <- matrix(runif(100,0,1),10,10)
colnames(net) <- rownames(net) <- letters[1:10]
node_level_covariates <- data.frame(Age = c(25,30,34,27,36,39,27,28,35,40),
	                                Height = c(70,70,67,58,65,67,64,74,76,80),
	                                Type = c("A","B","B","A","A","A","B","B","C","C"))
rownames(node_level_covariates) <- letters[1:10]
network_covariate <- net + matrix(rnorm(100,0,.5),10,10)
formula <- net ~ recip + edges + sender("Age") + 
netcov("network_covariate") + nodefactor("Type",base = "A")

test <- gergm(formula,
	          covariate_data = node_level_covariates,
	          network_is_directed = TRUE,
	          use_MPLE_only = FALSE,
	          estimation_method = "Metropolis",
	          maximum_number_of_lambda_updates = 5,
	          maximum_number_of_theta_updates = 5,
	          number_of_networks_to_simulate = 100000,
	          thin = 1/10,
	          proposal_variance = 0.5,
	          downweight_statistics_together = TRUE,
	          MCMC_burnin = 50000,
	          seed = 456,
	          convergence_tolerance = 0.01,
	          MPLE_gain_factor = 0,
	          force_x_theta_update = 2)

# Generate Estimate Plot
Estimate_Plot(test)
# Generate GOF Plot
GOF(test)
# Generate Trace Plot
Trace_Plot(test)
}
}