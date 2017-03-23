# Simulate a gergm
Simulate_GERGM <- function(GERGM_Object,
                           coef = GERGM_Object@theta.par,
                           seed1,
						               possible.stats,
						               verbose = TRUE,
						               parallel = FALSE,
						               predict_conditional_edges = FALSE,
						               i = NULL,
						               j = NULL) {
  # object: an object of class "gergm"

  sample_every <- floor(1/GERGM_Object@thin)
  thetas <- GERGM_Object@theta.par
  num.nodes <- GERGM_Object@num_nodes
  triples <- GERGM_Object@statistic_auxiliary_data$triples
  pairs <- GERGM_Object@statistic_auxiliary_data$pairs

  # if we are dealing with an undirected network
  undirect_network <- 0
  if (!GERGM_Object@directed_network) {
    undirect_network <- 1
  }

  # if we are dealing with a correlation network
  is_correlation_network <- 0
  if (GERGM_Object@is_correlation_network) {
    is_correlation_network <- 1
    undirect_network <- 1
  }

  # Gibbs Simulation
  if (GERGM_Object@estimation_method == "Gibbs") {
    nets <- Gibbs_Sampler(GERGM_Object,
                          thetas,
                          MCMC.burnin = GERGM_Object@burnin,
                          num.draws = GERGM_Object@number_of_simulations,
                          thin = GERGM_Object@thin,
                          start = NULL,
                          num.nodes = num.nodes,
                          directed = TRUE,
                          possible.stats = possible.stats)
    # Calculate the network statistics over all of the simulated networks
    for(i in 1:dim(nets)[3]) {
      temp <- calculate_h_statistics(
        GERGM_Object,
        GERGM_Object@statistic_auxiliary_data,
        all_weights_are_one = FALSE,
        calculate_all_statistics = TRUE,
        use_constrained_network = TRUE,
        network = nets[,,i])
      if (i == 1) {
        h.statistics <- matrix(0, nrow = dim(nets)[3], ncol = length(temp))
        h.statistics[i,] <- temp
      } else {
        h.statistics[i,] <- temp
      }
    }

    acceptance.rate <- NULL

    h.statistics <- as.data.frame(h.statistics)
    colnames(h.statistics) <- GERGM_Object@full_theta_names
  }

  # Metropolis Hastings Simulation
  if (GERGM_Object@estimation_method == "Metropolis") {

    # prepare variables for use with MH sampler
    store <- ceiling((GERGM_Object@number_of_simulations + GERGM_Object@burnin)/sample_every)
    nsim <- GERGM_Object@number_of_simulations + GERGM_Object@burnin
    dw <- as.numeric(GERGM_Object@downweight_statistics_together)

    # get the statistic auxiliary data list object
    sad <- GERGM_Object@statistic_auxiliary_data
    num_non_base_statistics <- sum(GERGM_Object@non_base_statistic_indicator)

    # prepare stochastic MH input
    num_unique_random_triad_samples <- 100
    smh <- generate_stochastic_MH_triples_pairs(
      GERGM_Object@stochastic_MH_proportion,
      GERGM_Object@use_stochastic_MH,
      triples = triples,
      pairs = pairs,
      samples = num_unique_random_triad_samples)

    # extract relevant quantites
    random_triad_samples <- smh$random_triples
    random_dyad_samples <- smh$random_pairs

    # do not use the multiplicative factor unless we are using stochastic MH
    if (GERGM_Object@use_stochastic_MH) {
      p_ratio_multaplicative_factor <- 1 / GERGM_Object@stochastic_MH_proportion
    } else {
      p_ratio_multaplicative_factor <- 1
    }

    rows_to_use <- sad$specified_rows_to_use - 1
    for (k in 1:length(rows_to_use)) {
      rows_to_use[k] <- max(rows_to_use[k], 0)
    }

    # if we are doing conditional edge prediction, then we only want to simulate
    # networks with one edge changing, conditional on the rest of the network
    # which will be fixed.
    if (predict_conditional_edges) {
      samples <- Individual_Edge_Conditional_Prediction(
        number_of_iterations = nsim,
        shape_parameter = GERGM_Object@proposal_variance,
        number_of_nodes = num.nodes,
        statistics_to_use = GERGM_Object@stats_to_use - 1,
        initial_network = GERGM_Object@bounded.network,
        take_sample_every = sample_every,
        thetas = thetas,
        triples = triples - 1,
        pairs = pairs - 1,
        alphas = GERGM_Object@weights,
        together = dw,
        seed = seed1,
        number_of_samples_to_store = store,
        using_correlation_network = is_correlation_network,
        undirect_network = undirect_network,
        parallel = parallel,
        use_selected_rows = sad$specified_selected_rows_matrix - 1,
        save_statistics_selected_rows_matrix = sad$full_selected_rows_matrix - 1,
        rows_to_use = rows_to_use,
        base_statistics_to_save = sad$full_base_statistics_to_save - 1,
        base_statistic_alphas = sad$full_base_statistic_alphas,
        num_non_base_statistics = num_non_base_statistics,
        non_base_statistic_indicator = GERGM_Object@non_base_statistic_indicator,
        p_ratio_multaplicative_factor = p_ratio_multaplicative_factor,
        random_triad_sample_list = random_triad_samples,
        random_dyad_sample_list = random_dyad_samples,
        use_triad_sampling = GERGM_Object@use_stochastic_MH,
        num_unique_random_triad_samples = num_unique_random_triad_samples,
        i = i - 1,
        j = j - 1)
    } else {
      # if we are not using the distribtuion estimator
      if (GERGM_Object@distribution_estimator == "none") {
        # take samples using MH
        samples <- Extended_Metropolis_Hastings_Sampler(number_of_iterations = nsim,
          shape_parameter = GERGM_Object@proposal_variance,
          number_of_nodes = num.nodes,
          statistics_to_use = GERGM_Object@stats_to_use - 1,
          initial_network = GERGM_Object@bounded.network,
          take_sample_every = sample_every,
          thetas = thetas,
          triples = triples - 1,
          pairs = pairs - 1,
          alphas = GERGM_Object@weights,
          together = dw,
          seed = seed1,
          number_of_samples_to_store = store,
          using_correlation_network = is_correlation_network,
          undirect_network = undirect_network,
          parallel = parallel,
          use_selected_rows = sad$specified_selected_rows_matrix - 1,
          save_statistics_selected_rows_matrix = sad$full_selected_rows_matrix - 1,
          rows_to_use = rows_to_use,
          base_statistics_to_save = sad$full_base_statistics_to_save - 1,
          base_statistic_alphas = sad$full_base_statistic_alphas,
          num_non_base_statistics = num_non_base_statistics,
          non_base_statistic_indicator = GERGM_Object@non_base_statistic_indicator,
          p_ratio_multaplicative_factor = p_ratio_multaplicative_factor,
          random_triad_sample_list = random_triad_samples,
          random_dyad_sample_list = random_dyad_samples,
          use_triad_sampling = GERGM_Object@use_stochastic_MH,
          num_unique_random_triad_samples = num_unique_random_triad_samples,
          include_diagonal = GERGM_Object@include_diagonal)
      } else {
        # if we are using the distribution estimator
        # take samples using MH

        # if FALSE, then we use the joint distribtuion sampler
        rowwise_distribution <- FALSE
        if (GERGM_Object@distribution_estimator == "rowwise-marginal") {
          rowwise_distribution <- TRUE
        }
        samples <- Distribution_Metropolis_Hastings_Sampler(
          number_of_iterations = nsim,
          variance = GERGM_Object@proposal_variance,
          number_of_nodes = num.nodes,
          statistics_to_use = GERGM_Object@stats_to_use - 1,
          initial_network = GERGM_Object@bounded.network,
          take_sample_every = sample_every,
          thetas = thetas,
          triples = triples - 1,
          pairs = pairs - 1,
          alphas = GERGM_Object@weights,
          together = dw,
          seed = seed1,
          number_of_samples_to_store = store,
          parallel = parallel,
          use_selected_rows = sad$specified_selected_rows_matrix - 1,
          save_statistics_selected_rows_matrix = sad$full_selected_rows_matrix - 1,
          rows_to_use = rows_to_use,
          base_statistics_to_save = sad$full_base_statistics_to_save - 1,
          base_statistic_alphas = sad$full_base_statistic_alphas,
          num_non_base_statistics = num_non_base_statistics,
          non_base_statistic_indicator = GERGM_Object@non_base_statistic_indicator,
          p_ratio_multaplicative_factor = p_ratio_multaplicative_factor,
          random_triad_sample_list = random_triad_samples,
          random_dyad_sample_list = random_dyad_samples,
          use_triad_sampling = GERGM_Object@use_stochastic_MH,
          num_unique_random_triad_samples = num_unique_random_triad_samples,
          rowwise_distribution = rowwise_distribution)
      }

    }

    # keep only the networks after the burnin
    start <- floor(GERGM_Object@burnin/sample_every) + 1
    end <- length(samples[[3]][,1])
    nets <- samples[[2]][, , start:end]
    # Note: these statistics will be the adjusted statistics (for use in the
    # MCMCMLE procedure)

    # more markov chain diagnostics
    average_log_prob_accept <- mean(samples[[5]])
    cat("Average log probability of accepting a proposal:",
        average_log_prob_accept,
        ".\nStandard deviation of log probability of accepting proposal:",
        sd(samples[[5]]),"\n")
    if (!is.finite(average_log_prob_accept)) {
      warning("It appears there is a problem with Metropolis Hastings, consider increasing proposal variance.")
    }

    GERGM_Object <- store_console_output(GERGM_Object,
      paste("Average log probability of accepting a proposal:",
            average_log_prob_accept,
            ".\nStandard deviation of log probability of accepting proposal:",
            sd(samples[[5]]),"\n"))

    average_edge_weight <- mean(samples[[4]])
    cat("Average (constrained) simulated network density:",
        average_edge_weight, "\n")
    GERGM_Object <- store_console_output(GERGM_Object,
      paste("Average (constrained) simulated network density:",
            average_edge_weight, "\n"))


    h.statistics <- samples[[3]][start:end,]
    acceptance.rate <- mean(samples[[1]])
    if (verbose) {
      cat("Metropolis Hastings Acceptance Rate (target = ",
          GERGM_Object@target_accept_rate," ): ",
          acceptance.rate, "\n", sep = "")
    }
    GERGM_Object <- store_console_output(GERGM_Object,
      paste("Metropolis Hastings Acceptance Rate (target = ",
            GERGM_Object@target_accept_rate,"):",
            acceptance.rate, "\n", sep = ""))

    P_Ratios = samples[[6]]
    Q_Ratios = samples[[7]]
    Proposed_Density = samples[[8]]
    Current_Density = samples[[9]]
    if (verbose) {
      cat("Average Q-Ratio:",mean(Q_Ratios),"Average P-Ratio:",mean(P_Ratios),
          "\nMean difference between proposed and current network densities:",
          mean(Proposed_Density - Current_Density ),"\n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,
     paste("Average Q-Ratio:",mean(Q_Ratios),"Average P-Ratio:",mean(P_Ratios),
           "\nMean difference between proposed and current network densities:",
           mean(Proposed_Density - Current_Density ),"\n", sep = ""))

    # make them a data frame and give them the correct row names
    h.statistics <- as.data.frame(h.statistics)
    colnames(h.statistics) <- GERGM_Object@full_theta_names
  }
  # if we are using MH, then return more diagnostics
  if (GERGM_Object@estimation_method == "Metropolis") {
    GERGM_Object@MCMC_output = list(Networks = nets,
                                    Statistics = h.statistics,
                                    Acceptance.rate = acceptance.rate,
                                    P_Ratios = samples[[6]],
                                    Q_Ratios = samples[[7]],
                                    Proposed_Density = samples[[8]],
                                    Current_Density = samples[[9]])
  } else {
    GERGM_Object@MCMC_output = list(Networks = nets,
                                    Statistics = h.statistics,
                                    Acceptance.rate = acceptance.rate)
  }
  return(GERGM_Object)
}
