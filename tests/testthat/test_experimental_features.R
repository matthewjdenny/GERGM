test_that("Hyperparameter optimization works", {
  skip_on_cran()
  skip("Skipping test that takes too long.")

  set.seed(12345)
  net <- matrix(runif(100,0,1),10,10)
  colnames(net) <- rownames(net) <- letters[1:10]
  node_level_covariates <- data.frame(Age = c(25,30,34,27,36,39,27,28,35,40),
                                      Height = c(70,70,67,58,65,67,64,74,76,80),
                                      Type = c("A","B","B","A","A","A","B","B","C","C"))
  rownames(node_level_covariates) <- letters[1:10]
  network_covariate <- net + matrix(rnorm(100,0,.5),10,10)

  formula <- net ~ edges + mutual + ttriads + sender("Age") +
    netcov("network_covariate") + nodematch("Type",base = "A")

  test <- gergm(formula,
                covariate_data = node_level_covariates,
                number_of_networks_to_simulate = 100000,
                thin = 1/100,
                proposal_variance = 0.1,
                MCMC_burnin = 50000,
                seed = 456,
                convergence_tolerance = 0.5,
                hyperparameter_optimization = TRUE)

  check_against <- c(1.319, -0.075, -0.017, -0.024,  3.091,  0.133, -1.837)
  check <- c(round(as.numeric(test@theta.coef[1,]),3),round(as.numeric(test@lambda.coef[1,]),3))
  expect_equal(check, check_against)

})





test_that("Model works for correlation networks", {
  skip_on_cran()

  set.seed(12345)
  #Function to generating a random positive-definite matrix with user-specified positive
  #eigenvalues
  # If eigenvalues are not specified, they are generated from a uniform
  #distribution

  Posdef <- function (n, ev = runif(n, 0, 10))
  {
    Z <- matrix(ncol=n, rnorm(n^2))
    decomp <- qr(Z)
    Q <- qr.Q(decomp)
    R <- qr.R(decomp)
    d <- diag(R)
    ph <- d / abs(d)
    O <- Q %*% diag(ph)
    Z <- t(O) %*% diag(ev) %*% O
    return(Z)
  }

  #Generate a random correlation matrix of dimension 10 x 10
  x <- rnorm(10)

  pdmat <- Posdef(n = 10)
  correlations <- pdmat / max(abs(pdmat))
  diag(correlations) <- 1
  net <- (correlations + t(correlations)) / 2
  colnames(net) <- rownames(net) <- letters[1:10]

  node_level_covariates <- data.frame(Type = c("A","B","B","A","A","A","B","B","C","C"))
  rownames(node_level_covariates) <- letters[1:10]


  formula <- net ~ edges + ttriads + nodematch("Type",base = "A")

  system.time({
    test <- gergm(formula,
                  covariate_data = node_level_covariates,
                  number_of_networks_to_simulate = 100000,
                  thin = 1/100,
                  proposal_variance = 0.2,
                  MCMC_burnin = 100000,
                  seed = 456,
                  convergence_tolerance = 0.5,
                  beta_correlation_model = TRUE,
                  convex_hull_proportion = 0.9,
                  convex_hull_convergence_proportion = 0.9,
                  sample_edges_at_a_time = 59,
                  use_previous_thetas = TRUE)
  })

  user  system elapsed
  38.320   0.624  43.158
  Theta Estimates:
    ttriads
  -0.02083058

  # make sure we are on observed scale:
  # test@MCMC_output$Networks[,,1]

  check_against <- c(-0.021,  0.070, -0.039, 40.285)
  check <- c(round(as.numeric(test@theta.coef[1,]),3),round(as.numeric(test@lambda.coef[1,]),3))
  expect_equal(check, check_against)

})


test_that("parallel threading works", {
  skip_on_cran()
  skip("We should not do threading on Travis")
  set.seed(12345)
  net <- matrix(runif(400,0,1),20,20)
  colnames(net) <- rownames(net) <- letters[1:20]
  node_level_covariates <- data.frame(Age = c(25,30,34,27,36,39,27,28,35,40,
                                              27,56,45,35,24,89,56,45,34,64),
                                      Height = c(70,70,67,58,65,67,64,74,76,80,
                                                 60,67,68,69,46,56,67,78,67,77),
                                      Type = c("A","B","B","A","A","A","B","B","C","C",
                                               "A","A","A","B","B","B","C","C","C","C"))
  rownames(node_level_covariates) <- letters[1:20]
  network_covariate <- net + matrix(rnorm(400,0,.5),20,20)

  formula <- net ~ edges + mutual + ttriads + sender("Age") +
    netcov("network_covariate") + nodematch("Type",base = "A")

  system.time({
    test <- gergm(formula,
                  covariate_data = node_level_covariates,
                  network_is_directed = TRUE,
                  use_MPLE_only = FALSE,
                  estimation_method = "Metropolis",
                  number_of_networks_to_simulate = 10000,
                  thin = 1/100,
                  proposal_variance = 0.1,
                  downweight_statistics_together = TRUE,
                  MCMC_burnin = 50000,
                  seed = 456,
                  convergence_tolerance = 0.01,
                  MPLE_gain_factor = 0,
                  force_x_theta_updates = 1,
                  hyperparameter_optimization = TRUE
    )
  })
  # user  system elapsed
  # 55.405   0.337  59.661

  system.time({
    test2 <- gergm(formula,
                  covariate_data = node_level_covariates,
                  network_is_directed = TRUE,
                  use_MPLE_only = FALSE,
                  estimation_method = "Metropolis",
                  number_of_networks_to_simulate = 10000,
                  thin = 1/100,
                  proposal_variance = 0.1,
                  downweight_statistics_together = TRUE,
                  MCMC_burnin = 50000,
                  seed = 456,
                  convergence_tolerance = 0.01,
                  MPLE_gain_factor = 0,
                  force_x_theta_updates = 1,
                  hyperparameter_optimization = TRUE,
                  parallel = TRUE,
                  cores = 2
    )
  })
  # user  system elapsed
  # 76.571   3.051  62.325



  check <- c(round(as.numeric(test@theta.coef[1,]),3),round(as.numeric(test@lambda.coef[1,]),3))
  check2 <- c(round(as.numeric(test2@theta.coef[1,]),3),round(as.numeric(test2@lambda.coef[1,]),3))
  expect_equal(check, check2)

})


test_that("sample at a time works", {
  skip_on_cran()
  set.seed(12345)
  net <- matrix(runif(400,0,1),20,20)
  colnames(net) <- rownames(net) <- letters[1:20]
  node_level_covariates <- data.frame(Age = c(25,30,34,27,36,39,27,28,35,40,
                                              27,56,45,35,24,89,56,45,34,64),
                                      Height = c(70,70,67,58,65,67,64,74,76,80,
                                                 60,67,68,69,46,56,67,78,67,77),
                                      Type = c("A","B","B","A","A","A","B","B","C","C",
                                               "A","A","A","B","B","B","C","C","C","C"))
  rownames(node_level_covariates) <- letters[1:20]
  network_covariate <- net + matrix(rnorm(400,0,.5),20,20)

  formula <- net ~ edges + mutual + ttriads + sender("Age") +
    netcov("network_covariate") + nodematch("Type",base = "A")

  system.time({
    test <- gergm(formula,
                  covariate_data = node_level_covariates,
                  network_is_directed = TRUE,
                  number_of_networks_to_simulate = 100000,
                  thin = 1/100,
                  proposal_variance = 0.1,
                  MCMC_burnin = 100000,
                  seed = 456,
                  convergence_tolerance = 0.5,
                  convex_hull_proportion = 0.9,
                  convex_hull_convergence_proportion = 0.9,
                  sample_edges_at_a_time = 229
    )
  })



  check <- c(round(as.numeric(test@theta.coef[1,]),3),round(as.numeric(test@lambda.coef[1,]),3))
  check2 <- c(round(as.numeric(test2@theta.coef[1,]),3),round(as.numeric(test2@lambda.coef[1,]),3))
  expect_equal(check, check2)

})


