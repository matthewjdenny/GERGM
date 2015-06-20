
# main function
gergm <- function(formula,
                  network_is_directed = c(TRUE, FALSE),
                  use_MPLE_only = c(FALSE, TRUE),
                  data_transformation = NULL,
                  estimation_method = c("Gibbs", "Metropolis"),
                  maximum_number_of_lambda_updates = 10,
                  maximum_number_of_theta_updates = 100,
                  number_of_networks_to_simulate = 500,
                  thin_chain_every_x = 1,
                  proposal_variance = 0.1,
                  exponential_weights = NULL,
                  downweight_statistics_together = TRUE,
                  MCMC_burnin = 100,
                  seed = 123,
                  convergence_tolerance = 0.01,
                  MPLE_gain_factor = 0.5){

  #' set logical values for whether we are using MPLE only, whether the network
  #' is directed, and which estimation method we are using
  use_MPLE_only <- use_MPLE_only[1] #default is FALSE
  network_is_directed <- network_is_directed[1] #default is TRUE
  estimation_method <- estimation_method[1] #default is Gibbs

  #' convert logical to numeric indicator
  if(downweight_statistics_together){
    downweight_statistics_together <- 1
  }else{
    downweight_statistics_together <- 0
  }

  #1. Create GERGM object from network


  #2. Estimate GERGM

  Model_Output <- Estimate_GERGM(formula,
                                 directed = network_is_directed,
                                 MPLE.only = use_MPLE_only,
                                 transform.data = data_transformation,
                                 method = estimation_method,
                                 max.num.iterations = maximum_number_of_lambda_updates,
                                 mc.num.iterations = maximum_number_of_theta_updates,
                                 nsim = number_of_networks_to_simulate,
                                 thin = thin_chain_every_x,
                                 shape.parameter = proposal_variance,
                                 exponential_weights = exponential_weights,
                                 together = downweight_statistics_together,
                                 MCMC.burnin = MCMC_burnin,
                                 seed = seed,
                                 tolerance = convergence_tolerance,
                                 gain.factor = MPLE_gain_factor)

  #3. Perform degeneracy diagnostics and create GOF plots


  #4. Return GERGM object


}
