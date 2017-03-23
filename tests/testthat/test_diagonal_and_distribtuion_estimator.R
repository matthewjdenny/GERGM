test_that("Diagonal model with no covariates runs", {
  skip_on_cran()
  skip("For time")
  ########################### 1. No Covariates #############################
  # Preparing an unbounded network without covariates for gergm estimation #
  #skip("Skipping test as it can only be run in the global environment.")

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

  check_against <- c(-2.016,  1.626,  0.118, -0.534)
  expect_equal(round(as.numeric(test@theta.coef[1,]),3), check_against)


})

test_that("Diagonal Model with covariates runs", {
  skip_on_cran()
  skip("For time")
  set.seed(12345)
  net <- matrix(runif(100,0,1),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]
  node_level_covariates <- data.frame(Age = c(25,30,34,27,36,39,27,28,35,40),
                                      Height = c(70,70,67,58,65,67,64,74,76,80),
                                      Type = c("A","B","B","A","A","A","B","B","C","C"))
  rownames(node_level_covariates) <- letters[1:10]
  network_covariate <- net + matrix(rnorm(100,0,.5),10,10)
  formula <- net ~ edges + mutual + ttriads + sender("Age") +
    netcov("network_covariate") + nodemix("Type",base = "A")

  test <- gergm(formula,
                covariate_data = node_level_covariates,
                network_is_directed = TRUE,
                use_MPLE_only = FALSE,
                estimation_method = "Metropolis",
                number_of_networks_to_simulate = 100000,
                thin = 1/100,
                proposal_variance = 0.5,
                downweight_statistics_together = TRUE,
                MCMC_burnin = 50000,
                seed = 456,
                convergence_tolerance = 0.01,
                MPLE_gain_factor = 0,
                force_x_theta_updates = 2)

  check_against <- c(0.824, -0.072, -0.016, -0.026, -0.024, -0.056, -0.055,
                     -0.035,  0.002, -0.040, -0.050,  3.061, 0.129, -1.931)
  expect_equal(c(round(as.numeric(test@theta.coef[1,]),3),round(as.numeric(test@lambda.coef[1,]),3)), check_against)


})


test_that("Row-wise distribution model runs", {
  skip_on_cran()
  skip("For time")
  ########################### 1. No Covariates #############################
  # Preparing an unbounded network without covariates for gergm estimation #
  #skip("Skipping test as it can only be run in the global environment.")

  set.seed(12345)
  net <- matrix(runif(100),10,10)
  for (i in 1:nrow(net)) {
    net[i,] <- net[i,]/sum(net[i,])
  }
  colnames(net) <- rownames(net) <- letters[1:10]

  # three parameter model
  formula <- net ~  edges(method = "endogenous") + mutual + ttriads

  test <- gergm(formula,
                normalization_type = "division",
                network_is_directed = TRUE,
                use_MPLE_only = FALSE,
                estimation_method = "Metropolis",
                number_of_networks_to_simulate = 200000,
                thin = 1/100,
                proposal_variance = 0.99,
                downweight_statistics_together = TRUE,
                MCMC_burnin = 100000,
                seed = 456,
                distribution_estimator = "rowwise-marginal",
                convergence_tolerance = 0.3)

  check_against <- c(-2.016,  1.626,  0.118, -0.534)
  expect_equal(round(as.numeric(test@theta.coef[1,]),3), check_against)


})
