run_mple <- function(GERGM_Object,
                     verbose,
                     seed2,
                     possible.stats,
                     outter_iteration_number = 1){

  # This function runs MPLE inside of the the MCMCMLE function

  statistics <- GERGM_Object@stats_to_use
  alphas <- GERGM_Object@weights
  together <- 1
  if (!GERGM_Object@downweight_statistics_together) {
    together <- 0
  }
  if (verbose) {
    cat("Estimating Initial Values for Theta via MPLE... \n")
  }
  GERGM_Object <- store_console_output(GERGM_Object,"Estimating Initial Values for Theta via MPLE... \n")

  if (GERGM_Object@is_correlation_network) {
    theta.init <- mple.corr(GERGM_Object@network,
                            GERGM_Object@bounded.network,
                            statistics = GERGM_Object@stats_to_use,
                            directed = GERGM_Object@directed_network,
                            alphas = rep(1, length(possible.stats)),
                            together = 1)
  } else if (GERGM_Object@weighted_MPLE) {
    # if we are using weighted MPLE with integration over edge weights then
    # allow exponential weighting in MPLE objective

    # if we are on the first iteration, use regular MPLE as an initialization.
    # Otherwise, use the theta estimates from the last round instead so we do
    # not have to worry about our thetas running off as much.

    # currently deprecating for now so that we can just work with weighted mple
    if (outter_iteration_number == 1) {
      theta.init <- ex_mple(GERGM_Object,
                            verbose = verbose)
    } else {
      theta.init <- list()
      theta.init$par <- GERGM_Object@theta.par
    }
    # update theta.init using weighted MPLE
    theta.init <- mple_weighted(GERGM_Object = GERGM_Object,
                       statistics = GERGM_Object@stats_to_use,
                       possible.stats = possible.stats,
                       verbose = verbose,
                       prev_ests = as.numeric(theta.init$par))

  } else {
    theta.init <- ex_mple(GERGM_Object, verbose = verbose)
  }
  if (verbose) {
    cat("\nMPLE Thetas: ", theta.init$par, "\n")
  }
  GERGM_Object <- store_console_output(GERGM_Object,
    paste("\nMPLE Thetas: ", theta.init$par, "\n"))
  num.nodes <- GERGM_Object@num_nodes
  triples <- t(combn(1:num.nodes, 3))
  if (GERGM_Object@is_correlation_network) {
    init.statistics <- calculate_h_statistics(
      GERGM_Object,
      GERGM_Object@statistic_auxiliary_data,
      all_weights_are_one = FALSE,
      calculate_all_statistics = TRUE,
      use_constrained_network = FALSE)
    obs.stats <- calculate_h_statistics(
      GERGM_Object,
      GERGM_Object@statistic_auxiliary_data,
      all_weights_are_one = FALSE,
      calculate_all_statistics = FALSE,
      use_constrained_network = FALSE)
  } else {
    init.statistics <- calculate_h_statistics(
      GERGM_Object,
      GERGM_Object@statistic_auxiliary_data,
      all_weights_are_one = FALSE,
      calculate_all_statistics = TRUE,
      use_constrained_network = TRUE)
    obs.stats <- calculate_h_statistics(
      GERGM_Object,
      GERGM_Object@statistic_auxiliary_data,
      all_weights_are_one = FALSE,
      calculate_all_statistics = FALSE,
      use_constrained_network = TRUE)
  }


  #cat("Observed Values of Selected Statistics:", "\n", obs.stats, "\n")
  ####################################################################
  GERGM_Object@reduced_weights <- alphas
  GERGM_Object@theta.par <- as.numeric(theta.init$par)

  # if we are not doing a fisher update
  theta <- list()
  theta$par <- as.numeric(theta.init$par)

  # if we are going to do a fisher update to MPLE thetas
  if (GERGM_Object@MPLE_gain_factor > 0) {
    cat("Updating MPLE estimates via a one step Fisher update...\n")
    GERGM_Object <- Simulate_GERGM(GERGM_Object,
                                   seed1 = seed2,
                                   possible.stats = possible.stats,
                                   verbose = verbose)
    indicies <- GERGM_Object@statistic_auxiliary_data$specified_statistic_indexes_in_full_statistics

    hsn <- GERGM_Object@MCMC_output$Statistics[,indicies]

    #Calculate covariance estimate (to scale initial guess theta.init)
    z.bar <- NULL
    if (class(hsn) == "numeric") {
      hsn <- matrix(hsn,ncol = 1,nrow = length(hsn))
      z.bar <- sum(hsn)
    } else {
      z.bar <- colSums(hsn)
    }

    #cat("z.bar", "\n", z.bar, "\n")
    Cov.est <- 0
    for (i in 1:dim(hsn)[1]) {
      Cov.est <- matrix(as.numeric(hsn[i,]), ncol = 1) %*%
        t(matrix(as.numeric(hsn[i,]), ncol = 1)) + Cov.est
    }
    Cov.est <- (Cov.est) - z.bar %*% t(z.bar)
    #cat("Cov.est", "\n", Cov.est)

    # try to update but if we get a zero percent accept rate then
    # do not udpate.
    try({
      D.inv <- solve(Cov.est)
      #calculate
      theta$par <- as.numeric(theta.init$par) - GERGM_Object@MPLE_gain_factor *
        D.inv %*% (z.bar - obs.stats)
    })
    if (verbose) {
      cat("Adjusted Initial Thetas After Fisher Update:",theta$par, "\n\n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,paste("Adjusted Initial Thetas After Fisher Update:",theta$par, "\n\n"))
  }


  # if we are using the experimental optimization feature.
  if (GERGM_Object@hyperparameter_optimization) {
    if (GERGM_Object@using_grid_optimization) {
      cat("Optimizing thetas, this may take days...\n")
      theta$par <- optimize_initialization(GERGM_Object,
                                           verbose,
                                           seed2,
                                           possible.stats,
                                           theta,
                                           statistics)
      cat("Updated thetas after grid search:",theta$par,"\n\n")
    }
  }

  theta$par <- as.numeric(theta$par)
  return(list(GERGM_Object = GERGM_Object,
              theta = theta,
              statistics = statistics,
              init.statistics = init.statistics,
              hessian = theta.init$hessian))
}
