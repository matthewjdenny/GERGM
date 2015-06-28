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
########################### 1. No Covariates #############################
# Preparing an unbounded network without covariates for gergm estimation #
net <- matrix(rnorm(100,0,20),10,10)
colnames(net) <- rownames(net) <- letters[1:10]
network <- Prepare_Network_and_Covariates(raw_network = net,
                                          normalization_type = "division")

formula <- "network ~ recip + edgeweight"

test <- gergm(formula,
              network_is_directed = TRUE,
              use_MPLE_only = FALSE,
              data_transformation = NULL,
              estimation_method = "Gibbs",
              maximum_number_of_lambda_updates = 1,
              maximum_number_of_theta_updates = 5,
              number_of_networks_to_simulate = 40000,
              thin = 1/10,
              proposal_variance = 0.5,
              exponential_weights = NULL,
              downweight_statistics_together = TRUE,
              MCMC_burnin = 10000,
              seed = 456,
              convergence_tolerance = 0.01,
              MPLE_gain_factor = 0,
              force_x_theta_update = 4,
              output_directory = getwd(),
              output_name= "Testing")

########################### 2. Covariates #############################
# Preparing an unbounded network with covariates for gergm estimation #
net <- matrix(rnorm(100,0,20),10,10)
# necessary to identify nodes
colnames(net) <- rownames(net) <- letters[1:10]
node_level_covariates <- data.frame(Age = c(25,30,34,27,36,39,27,28,35,40),
                                    Height = c(70,70,67,58,65,67,64,74,76,80))
# necessary to identify nodes									
rownames(node_level_covariates) <- letters[1:10]
network_covariate <- net + matrix(rnorm(100,0,.5),10,10) %*% matrix(1:100,10,10)
output <- Prepare_Network_and_Covariates(raw_network = net,
                                         covariate_data = node_level_covariates,
                                         type_of_effect = c("sender","both"),
                                         network_covariates = network_covariate,
                                         network_covariate_names = "My_Net")

transformed_covariates <- output$transformed_covariates
network <- output$network

formula <- "network ~ recip + edgeweight"

#note that this model is currently degenerate -- illustrating automatic stopping
test <- gergm(formula,
              network_is_directed = TRUE,
              use_MPLE_only = FALSE,
              data_transformation = transformed_covariates,
              estimation_method = "Metropolis",
              maximum_number_of_lambda_updates = 1,
              maximum_number_of_theta_updates = 5,
              number_of_networks_to_simulate = 100000,
              thin = 1/10,
              proposal_variance = 0.5,
              exponential_weights = NULL,
              downweight_statistics_together = TRUE,
              MCMC_burnin = 10000,
              seed = 456,
              convergence_tolerance = 0.01,
              MPLE_gain_factor = 0,
              force_x_theta_update = 1,
              output_directory = getwd(),
              output_name= "Testing")
}