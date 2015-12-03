test_that("Simple model with no covariates runs", {
  skip_on_cran()
  ########################### 1. No Covariates #############################
  # Preparing an unbounded network without covariates for gergm estimation #
  set.seed(12345)
  net <- matrix(rnorm(100,0,20),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]
  formula <- net ~ recip + edges + ttriads

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
                force_x_theta_updates = 4)

formula2 <- net ~ recip + edges + ttriads + in2star + ctriads

test <- gergm(formula2,
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
              force_x_theta_updates = 4)
})

test_that("Model with covariates runs", {
  skip_on_cran()
  set.seed(12345)
  net <- matrix(runif(100,0,1),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]
  node_level_covariates <- data.frame(Age = c(25,30,34,27,36,39,27,28,35,40),
                                      Height = c(70,70,67,58,65,67,64,74,76,80),
                                      Type = c("A","B","B","A","A","A","B","B","C","C"))
  rownames(node_level_covariates) <- letters[1:10]
  network_covariate <- net + matrix(rnorm(100,0,.5),10,10)
  formula <- net ~ recip + edges + ttriads + sender("Age") +
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
                force_x_theta_updates = 2)
})
