run_mple <- function(GERGM_Object,
                     verbose,
                     seed2,
                     possible.stats){

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
    theta.init <- mple.corr(GERGM_Object@network, GERGM_Object@bounded.network,
                            statistics = GERGM_Object@stats_to_use,
                            directed = GERGM_Object@directed_network,
                            alphas = rep(1, length(possible.stats)),
                            together = 1)
  } else if (GERGM_Object@weighted_MPLE) {
    # if we are using weighted MPLE with integration over edge weights then
    # allow exponential weighting in MPLE objective
    theta.init <- mple(GERGM_Object@bounded.network,
                       statistics = GERGM_Object@stats_to_use,
                       directed = GERGM_Object@directed_network,
                       alphas = alphas,
                       together = together,
                       weighted_MPLE = TRUE,
                       verbose = verbose)

  } else {
    theta.init <- mple(GERGM_Object@bounded.network,
                       statistics = GERGM_Object@stats_to_use,
                       directed = GERGM_Object@directed_network,
                       alphas = rep(1, length(possible.stats)),
                       together = 1,
                       weighted_MPLE = FALSE,
                       verbose = verbose)
  }
  if (verbose) {
    cat("\nMPLE Thetas: ", theta.init$par, "\n")
  }
  GERGM_Object <- store_console_output(GERGM_Object, paste("\nMPLE Thetas: ", theta.init$par, "\n"))
  num.nodes <- GERGM_Object@num_nodes
  triples <- t(combn(1:num.nodes, 3))
  if (GERGM_Object@is_correlation_network) {
    # initialize the network with the observed network
    initial_network <- GERGM_Object@network
    # calculate the statistics of the original network
    init.statistics <- h2(GERGM_Object@network,
                          triples = triples,
                          statistics = rep(1, length(possible.stats)),
                          alphas = alphas,
                          together = GERGM_Object@downweight_statistics_together)
    obs.stats <- h2(GERGM_Object@network,
                    triples = triples,
                    statistics = GERGM_Object@stats_to_use,
                    alphas = alphas,
                    together = GERGM_Object@downweight_statistics_together)
  } else {
    # initialize the network with the observed network
    initial_network <- GERGM_Object@bounded.network
    # calculate the statistics of the original network
    init.statistics <- h2(GERGM_Object@bounded.network,
                          triples = triples,
                          statistics = rep(1, length(possible.stats)),
                          alphas = alphas,
                          together = GERGM_Object@downweight_statistics_together)
    obs.stats <- h2(GERGM_Object@bounded.network,
                    triples = triples,
                    statistics = GERGM_Object@stats_to_use,
                    alphas = alphas,
                    together = GERGM_Object@downweight_statistics_together)
  }


  #cat("Observed Values of Selected Statistics:", "\n", obs.stats, "\n")
  ####################################################################
  alps <- alphas[which(statistics == 1)]
  GERGM_Object@reduced_weights <- alps
  GERGM_Object@theta.par <- theta.init$par

  # if we are not doing a fisher update
  theta <- theta.init

  # if we are going to do a fisher update to MPLE thetas
  if (GERGM_Object@MPLE_gain_factor > 0) {
    GERGM_Object <- Simulate_GERGM(GERGM_Object,
                                   seed1 = seed2,
                                   possible.stats = possible.stats,
                                   verbose = verbose)

    hsn <- GERGM_Object@MCMC_output$Statistics[,which(GERGM_Object@stats_to_use == 1)]

    #Calculate covariance estimate (to scale initial guess theta.init)
    z.bar <- NULL
    if (class(hsn) == "numeric") {
      hsn <- matrix(hsn,ncol = 1,nrow = length(hsn))
      z.bar <- sum(hsn) / 20
    } else {
      z.bar <- colSums(hsn) / 20
    }

    #cat("z.bar", "\n", z.bar, "\n")
    Cov.est <- 0
    for (i in 1:dim(hsn)[1]) {
      Cov.est <- matrix(as.numeric(hsn[i,]), ncol = 1) %*%
        t(matrix(as.numeric(hsn[i,]), ncol = 1)) + Cov.est
    }
    Cov.est <- (Cov.est / 20) - z.bar %*% t(z.bar)
    #cat("Cov.est", "\n", Cov.est)

    # try to update but if we get a zero percent accept rate then
    # do not udpate.
    try({
      D.inv <- solve(Cov.est)
      #calculate
      theta$par <- theta.init$par - GERGM_Object@MPLE_gain_factor *
        D.inv %*% (z.bar - obs.stats)
    })
    if (verbose) {
      cat("Adjusted Initial Thetas After Fisher Update:",theta$par, "\n\n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,paste("Adjusted Initial Thetas After Fisher Update:",theta$par, "\n\n"))
  }

  return(list(GERGM_Object = GERGM_Object,
              theta = theta,
              statistics = statistics,
              init.statistics = init.statistics))
}
