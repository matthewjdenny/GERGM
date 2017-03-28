test_that("calculating statistics based on a categorical node level covariate works", {
  skip_on_cran()
  skip("For time")

  set.seed(12345)
  net <- matrix(runif(100,0,1),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]
  node_level_covariates <- data.frame(Age = c(25,30,34,27,36,39,27,28,35,40),
                                      Height = c(70,70,67,58,65,67,64,74,76,80),
                                      Type = c("A","B","B","A","A","C","B","B","C","C"))
  rownames(node_level_covariates) <- letters[1:10]
  network_covariate <- net + matrix(rnorm(100,0,.5),10,10)

  formula <- net ~ edges + mutual + ctriads + out2stars(covariate = "Type") + in2stars(covariate = "Type", base = "C") + sender("Age") +
    netcov("network_covariate") + nodematch("Type",base = "A")

  test <- gergm(formula,
                covariate_data = node_level_covariates,
                network_is_directed = TRUE,
                use_MPLE_only = FALSE,
                estimation_method = "Metropolis",
                number_of_networks_to_simulate = 100000,
                thin = 1/100,
                proposal_variance = 0.1,
                downweight_statistics_together = TRUE,
                MCMC_burnin = 50000,
                seed = 456,
                convergence_tolerance = 0.5,
                hyperparameter_optimization = TRUE
  )

  check_against <- c(0.706,  0.372, -3.749, -0.113, -6.246,  3.601, -3.516,
                     -0.013,  0.085,  2.501,  0.180, -2.026)
  check <- c(round(as.numeric(test@theta.coef[1,]),3),
             round(as.numeric(test@lambda.coef[1,]),3))
  expect_equal(check, check_against)

})




test_that("stochastic aproximation works", {
  skip_on_cran()
  skip("Still in development")
  # have not messed with this in a while (3-27-17) may not work.

  # try some simulations
  system.time({
    set.seed(12345)
    net <- matrix(runif(40000),200,200)
    diag(net) <- 0
    colnames(net) <- rownames(net) <- 1:200
    formula <- net ~ edges + ttriads + in2stars

    test <- simulate_networks(formula,
                              thetas = c(-.4,.3),
                              network_is_directed = TRUE,
                              simulation_method = "Metropolis",
                              number_of_networks_to_simulate = 100,
                              thin = 1/10,
                              proposal_variance = 0.0001,
                              downweight_statistics_together = TRUE,
                              MCMC_burnin = 50,
                              omit_intercept_term = TRUE,
                              seed = 456)
  })

  ncol(combn(1:500,3))
  system.time({
    set.seed(12345)
    net <- matrix(runif(250000),500,500)
    diag(net) <- 0
    colnames(net) <- rownames(net) <- 1:500
    formula <- net ~ edges + ttriads + in2stars

    test2 <- simulate_networks(formula,
                              thetas = c(-.4,.3),
                              network_is_directed = TRUE,
                              simulation_method = "Metropolis",
                              number_of_networks_to_simulate = 100,
                              thin = 1/10,
                              proposal_variance = 0.00005,
                              downweight_statistics_together = TRUE,
                              MCMC_burnin = 50,
                              omit_intercept_term = TRUE,
                              seed = 456,
                              use_stochastic_MH = TRUE,
                              stochastic_MH_proportion = 0.00004)
  })



  set.seed(12345)
  net <- matrix(runif(2500,0,1),50,50)
  colnames(net) <- rownames(net) <- 1:50

  formula <- net ~ edges + mutual + ttriads

  test <- gergm(formula,
                estimation_method = "Metropolis",
                number_of_networks_to_simulate = 1000,
                thin = 1/10,
                proposal_variance = 0.1,
                downweight_statistics_together = TRUE,
                MCMC_burnin = 500,
                seed = 456,
                convergence_tolerance = 0.5,
                hyperparameter_optimization = TRUE,
                use_stochastic_MH = TRUE,
                stochastic_MH_proportion = 0.1
  )

  check_against <- c(0.657,  0.212, -2.134,  0.234, -3.997,  2.960, -2.605,
                     -0.013,  0.057,  2.586, 0.170, -1.944)
  check <- c(round(as.numeric(test@theta.coef[1,]),3),round(as.numeric(test@lambda.coef[1,]),3))
  expect_equal(check, check_against)

})
