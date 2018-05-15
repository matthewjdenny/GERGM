test_that("Simple model with no covariates runs", {
  skip_on_cran()
  skip("For time")
  ########################### 1. No Covariates #############################
  # Preparing an unbounded network without covariates for gergm estimation #
  #skip("Skipping test as it can only be run in the global environment.")

  set.seed(12345)
  net <- matrix(rnorm(100,0,20),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]

  # one parameter model
  formula <- net ~ edges + ttriads

  test <- gergm(formula,
                number_of_networks_to_simulate = 40000,
                thin = 1/40,
                proposal_variance = 0.5,
                MCMC_burnin = 10000,
                seed = 456,
                convergence_tolerance = 0.5)

  check_against <- c(-0.081)
  expect_equal(round(as.numeric(test@theta.coef[1,]),3), check_against)
})
test_that("3 param model with no covariates runs", {
  skip_on_cran()
  skip("Weird travis errors")

  set.seed(12345)
  net <- matrix(rnorm(100,0,20),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]

  # three parameter model
  formula <- net ~  edges + mutual + ttriads

  test <- gergm(formula,
                number_of_networks_to_simulate = 40000,
                thin = 1/40,
                proposal_variance = 0.5,
                MCMC_burnin = 10000,
                seed = 456,
                convergence_tolerance = 0.5)

  check_against <- c(2.180, -0.268)
  expect_equal(round(as.numeric(test@theta.coef[1,]),3), check_against)

})
test_that("4 param model with no covariates runs", {
  skip_on_cran()
  skip("For time")

  set.seed(12345)
  net <- matrix(rnorm(100,0,20),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]

  # five parameter model
  formula2 <- net ~  edges + mutual + ttriads + in2stars

  test <- gergm(formula2,
              number_of_networks_to_simulate = 40000,
              thin = 1/40,
              proposal_variance = 0.5,
              downweight_statistics_together = TRUE,
              MCMC_burnin = 10000,
              seed = 456,
              convergence_tolerance = 0.5)

check_against <- c(2.329,  0.064, -0.498)
expect_equal(round(as.numeric(test@theta.coef[1,]),3), check_against)

})
test_that("Unidrected model with no covariates runs", {
  skip_on_cran()
  skip("Weird travis errors")
#check that code works for undirected network

  set.seed(12345)
  net <- matrix(rnorm(100,0,20),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]

  formula <- net ~ edges + ttriads + twostars

  test <- gergm(formula,
              network_is_directed = FALSE,
              number_of_networks_to_simulate = 40000,
              thin = 1/40,
              proposal_variance = 0.5,
              MCMC_burnin = 10000,
              seed = 456,
              convergence_tolerance = 0.5)

  check_against <- c(0.120, -0.282)
  expect_equal(round(as.numeric(test@theta.coef[1,]),3), check_against)

})
test_that("MPLE Only", {
  skip_on_cran()
#check that code works with MPLE only

  set.seed(12345)
  net <- matrix(rnorm(100,0,20),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]

  formula <- net ~ edges + ttriads + in2stars

  test <- gergm(formula,
                use_MPLE_only = TRUE,
                estimation_method = "Metropolis",
                number_of_networks_to_simulate = 40000,
                thin = 1/40,
                proposal_variance = 0.5,
                downweight_statistics_together = TRUE,
                MCMC_burnin = 10000,
                seed = 456,
                convergence_tolerance = .5)

  check_against <- c(0.187, -0.451)
  expect_equal(round(as.numeric(test@theta.coef[1,]),3), check_against)

})

test_that("Model with covariates runs", {
  skip_on_cran()
  skip("Weird travis errors")

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
                number_of_networks_to_simulate = 100000,
                thin = 1/100,
                proposal_variance = 0.5,
                MCMC_burnin = 50000,
                seed = 456,
                convergence_tolerance = 0.5)

  check_against <- c(0.764, -0.071, -0.016, -0.025, -0.023, -0.056, -0.056, -0.034,
                     0.002, -0.039, -0.050,  3.096,  0.128, -1.933)
  expect_equal(c(round(as.numeric(test@theta.coef[1,]),3),round(as.numeric(test@lambda.coef[1,]),3)), check_against)

})

test_that("Additional Model with covariates runs", {
  skip_on_cran()
  skip("Time")
  set.seed(12345)
  net <- matrix(runif(100,0,1),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]
  node_level_covariates <- data.frame(Age = c(25,30,34,27,36,39,27,28,35,40),
                                      Height = c(70,70,67,58,65,67,64,74,76,80),
                                      Type = c("A","B","B","A","A","A","B","B","C","C"))
  rownames(node_level_covariates) <- letters[1:10]
  network_covariate <- net + matrix(rnorm(100,0,.5),10,10)

  formula <- net ~ edges + mutual + ttriads + sender("Age") +
    netcov("network_covariate") + nodematch("Type")

  test <- gergm(formula,
                covariate_data = node_level_covariates,
                number_of_networks_to_simulate = 100000,
                thin = 1/100,
                proposal_variance = 0.5,
                MCMC_burnin = 50000,
                seed = 456,
                convergence_tolerance = 0.5,
                convex_hull_proportion = 0.9)

  check_against <- c(1.339, -0.074, -0.017, -0.024,  3.098,  0.132, -1.837)
  expect_equal(c(round(as.numeric(test@theta.coef[1,]),3),round(as.numeric(test@lambda.coef[1,]),3)), check_against)

})
