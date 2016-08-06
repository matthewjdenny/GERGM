test_that("SlackR integration works", {
  skip_on_cran()
  skip("We should not be pushing to slack automatically")

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


  test_webhook_url = "https://hooks.slack.com/services/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"

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
                force_x_theta_updates = 2,
                slackr_integration_list = list(model_name = "descriptive model name",
                     channel = "#slackrtesting",
                     incoming_webhook_url = test_webhook_url))

})





