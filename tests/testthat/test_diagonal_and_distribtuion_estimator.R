  test_that("Diagonal model with no covariates runs", {
  skip_on_cran()
    skip("For time...")

  set.seed(12345)
  net <- matrix(runif(100),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]

  # three parameter model
  formula <- net ~  edges(method = "endogenous") + mutual + ttriads

  test <- gergm(formula,
                normalization_type = "division",
                network_is_directed = TRUE,
                use_MPLE_only = FALSE,
                estimation_method = "Metropolis",
                number_of_networks_to_simulate = 40000,
                thin = 1/40,
                proposal_variance = 0.1,
                downweight_statistics_together = TRUE,
                MCMC_burnin = 10000,
                seed = 456,
                include_diagonal = TRUE)

  check_against <- c(-2.102,  1.475,  0.144, -1.064)
  expect_equal(round(as.numeric(test@theta.coef[1,]),3), check_against)


})



test_that("Row-wise distribution model runs", {
  skip_on_cran()
  skip("For time")

  set.seed(12345)
  net <- matrix(runif(100),10,10)
  for (i in 1:nrow(net)) {
    net[i,] <- net[i,]/sum(net[i,])
  }
  colnames(net) <- rownames(net) <- letters[1:10]

  # three parameter model
  formula <- net ~ mutual + ttriads

  test <- gergm(formula,
                number_of_networks_to_simulate = 200000,
                thin = 1/100,
                proposal_variance = 0.3,
                MCMC_burnin = 100000,
                seed = 456,
                distribution_estimator = "rowwise-marginal",
                convergence_tolerance = 0.5,
                parallel = TRUE,
                cores = 2,
                hyperparameter_optimization = TRUE,
                integration_intervals = 150,
                distribution_mple_regularization_weight = 0.02,
                force_x_theta_updates = 3)

  check_against <- c(1.311, -13.889, 5.779)
  expect_equal(round(as.numeric(test@theta.coef[1,]),3), check_against)


})
