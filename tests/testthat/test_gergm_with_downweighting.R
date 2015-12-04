test_that("Simple model with no covariates runs", {
  skip_on_cran()
  ########################### 1. No Covariates #############################
  # Preparing an unbounded network without covariates for gergm estimation #
  skip("Skipping test as it can only be run in the global environment.")

  set.seed(12345)
  net <- matrix(rnorm(100,0,20),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]


  # three parameter model
  formula <- net ~  edges + mutual(0.95) +  ttriads(0.95)

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
                MPLE_gain_factor = 0.05,
                force_x_theta_updates = 4)

  check_against <- c(2.197, -0.023, -0.395)
  expect_equal(round(as.numeric(test@theta.coef[1,]),3), check_against)


  #check that code works for undirected network

  formula <- net ~ edges(alpha = 0.8) + ttriads(alpha = 0.9) + twostars(0.9)

  test <- gergm(formula,
                normalization_type = "division",
                network_is_directed = FALSE,
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
                MPLE_gain_factor = 0.08,
                force_x_theta_updates = 4)


  set.seed(12345)
  net <- matrix(runif(100,0,1),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]
  node_level_covariates <- data.frame(Age = c(25,30,34,27,36,39,27,28,35,40),
                                      Height = c(70,70,67,58,65,67,64,74,76,80),
                                      Type = c("A","B","B","A","A","A","B","B","C","C"))
  rownames(node_level_covariates) <- letters[1:10]
  network_covariate <- net + matrix(rnorm(100,0,.5),10,10)
  formula <- net ~ edges(0.95) + mutual(0.95) + ttriads(0.95) + sender("Age") +
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
                MPLE_gain_factor = 0.05,
                force_x_theta_updates = 2)

  check_against <- c(2.546, -0.101, -0.743,  0.583, -0.227, -0.043, -0.040,  0.105, -1.877)
  check <- c(round(as.numeric(test@theta.coef[1,]),3),round(as.numeric(test@lambda.coef[1,]),3))
  expect_equal(check, check_against)

})