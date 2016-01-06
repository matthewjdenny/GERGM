hysteresis <- function(GERGM_Object,
                       networks_to_simulate = 1000,
                       burnin = 500,
                       range = NULL,
                       steps = 50,
                       simulation_method = c("Gibbs", "Metropolis"),
                       proposal_variance = 0.1,
                       seed = 12345,
                       thin = 1,
                       output_directory = NULL,
                       output_name = NULL,
                       parallel = TRUE,
                       cores = 1){

  # preliminaries
  possible_structural_terms <- c("out2stars",
                                 "in2stars",
                                 "ctriads",
                                 "mutual",
                                 "ttriads",
                                 "edges")
  currentwd <- getwd()


  # figure out how many statistics we need to simulate values for
  num_network_terms <- length(GERGM_Object@theta.par)
  Hysteresis_Results <- vector(mode = "list", length = num_network_terms)

  for(i in 1:num_network_terms){
    # figure out the range of values for each parameter


    GERGM_Object <- Simulate_GERGM(
      GERGM_Object,
      nsim = networks_to_simulate,
      method = simulation_method,
      MCMC.burnin = burnin,
      thin = thin,
      shape.parameter = proposal_variance,
      together = GERGM_Object@downweight_statistics_together,
      seed1 = seed,
      possible.stats = possible_structural_terms)
  }


  # clean up a return everything
  setwd(currentwd)
  return(Hysteresis_Results)
}
