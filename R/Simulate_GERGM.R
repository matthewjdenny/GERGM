# Simulate a gergm
Simulate_GERGM <- function(GERGM_Object,
                           nsim,
                           method,
                           MCMC.burnin,
                           weights = GERGM_Object@weights,
                           coef = GERGM_Object@theta.par,
                           thin ,
                           shape.parameter ,
                           together ,
                           seed1,
						               possible.stats ) {
  # object: an object of class "gergm"

  sample_every <- floor(1/thin)
  thetas <- GERGM_Object@theta.par
  num.nodes <- GERGM_Object@num_nodes
  triples <- t(combn(1:num.nodes, 3))
  pairs <- t(combn(1:num.nodes, 2))

  # Gibbs Simulation
  if (method == "Gibbs") {
    nets <- Gibbs_Sampler(GERGM_Object,
                          thetas,
                          MCMC.burnin = MCMC.burnin,
                          num.draws = nsim,
                          thin = thin,
                          start = NULL,
                          num.nodes = num.nodes,
                          directed = TRUE,
                          possible.stats = possible.stats)
    # Calculate the network statistics over all of the simulated networks
    h.statistics <- t(apply(nets, 3, h2, triples = triples,
                            statistics = rep(1, 6), alphas = rep(1, 6),
                            together = together))
    acceptance.rate <- NULL
  }

  # Metropolis Hastings Simulation
  if (method == "Metropolis") {
    #need to put the thetas into a full length vector for MH function
    stat.indx <- which(GERGM_Object@stats_to_use > 0)
    #cat("stat.idx",stat.indx,"\n" )
    full_thetas <- rep(0, length(GERGM_Object@stats_to_use))
    for (i in 1:length(thetas)) {
      full_thetas[stat.indx[i]] <- thetas[i]
    }
    #cat("Current Theta Estimates:",thetas,"\n")
    store <- ceiling((nsim + MCMC.burnin)/sample_every)
    samples <- Metropolis_Hastings_Sampler(
      number_of_iterations = nsim + MCMC.burnin,
      shape_parameter = shape.parameter,
      number_of_nodes = num.nodes,
      statistics_to_use = GERGM_Object@stats_to_use,
      initial_network = GERGM_Object@bounded.network,
      take_sample_every = sample_every,
      thetas = full_thetas,
      triples = triples - 1,
      pairs = pairs - 1,
      alphas = GERGM_Object@weights,
      together = together,
      seed = seed1,
      number_of_samples_to_store = store)
    # keep only the networks after the burnin
    start <- floor(MCMC.burnin/sample_every) + 1
    end <- length(samples[[3]][,1])
    nets <- samples[[2]][, , start:end]
    # Note: these statistics will be the adjusted statistics (for use in the
    # MCMCMLE procedure)

    h.statistics <- samples[[3]][start:end,]
    acceptance.rate <- mean(samples[[1]])
    cat("Metropolis Hastings Acceptance Rate (target = 0.25):", acceptance.rate, "\n")
    GERGM_Object <- store_console_output(GERGM_Object,paste("Metropolis Hastings Acceptance Rate (target = 0.25):", acceptance.rate, "\n"))

  }
  h.statistics = data.frame(out2stars = h.statistics[, 1],
                            in2stars = h.statistics[, 2],
                            ctriads = h.statistics[, 3],
                            recip = h.statistics[, 4],
                            ttriads = h.statistics[, 5],
                            edgeweight = h.statistics[, 6])

  GERGM_Object@MCMC_output = list(Networks = nets,
                            Statistics = h.statistics,
                            Acceptance.rate = acceptance.rate)
  return(GERGM_Object)
}
