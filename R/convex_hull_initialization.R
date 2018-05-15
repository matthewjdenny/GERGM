convex_hull_initialization <- function(GERGM_Object,
                                       seed,
                                       possible.stats,
                                       verbose) {

  # start by calculating those statistics we are actually estimating for the
  # observed network. We will then compare these against the cloud of statistics
  # we simulate to see if the observed statistics live inside of the cloud.
  observed_network_stats <- calculate_h_statistics(
    GERGM_Object,
    GERGM_Object@statistic_auxiliary_data,
    all_weights_are_one = FALSE,
    calculate_all_statistics = FALSE,
    use_constrained_network = TRUE)

  # this is a logical indicating that we have not found an initialization where
  # the observed statistics are inside the convex hull yet.
  converged <- FALSE
  iteration <- 1

  cat("Starting thetas:",GERGM_Object@theta.par,"\n")

  while (!converged) {
    cat("Current iteration of convex hull initialization:",iteration,"\n")
    iteration <- iteration + 1

    # In each iteration, we first simulate a bunch of networks from the current
    # value for theta:
    GERGM_Object <- Simulate_GERGM(
      GERGM_Object,
      coef = GERGM_Object@theta.par,
      seed1 = seed,
      possible.stats = possible.stats,
      verbose = verbose,
      parallel = GERGM_Object@parallel_statistic_calculation)

    # now extract out the statistics for each simulated network
    sad <- GERGM_Object@statistic_auxiliary_data
    hsn <- GERGM_Object@MCMC_output$Statistics[,sad$specified_statistic_indexes_in_full_statistics]

    # deal with case where we only have one statistic, make sure we end up with
    # a matrix:
    if (class(hsn) == "numeric") {
      hsn <- matrix(hsn,
                    ncol = 1,
                    nrow = length(hsn))
    }

    # now we find the centroid of these statistics
    simulated_center <- colMeans(hsn)

    # now we find the two points closest to the observed network in the
    # statistic space:
    observed_distances <- rep(0,nrow(hsn))
    sim_center_distances <- rep(0,nrow(hsn))

    # basic function to calculate distances:
    calc_dists <- function(target,
                           test) {
      dist <- 0
      for (j in 1:length(target)) {
        dist <- dist + (target[j] - test[j])^2
      }
      dist <- sqrt(dist)
      return(dist)
    }

    # now loop through and calculate the distances:
    for (i in 1:nrow(hsn)) {
      observed_distances[i] <- calc_dists(observed_network_stats, hsn[i,])
      sim_center_distances[i] <- calc_dists(simulated_center, hsn[i,])
    }

    if (class(observed_distances) == "list") {
      observed_distances <- as.numeric(unlist(observed_distances))
    }
    # find the closest points
    ordering <- order(observed_distances, decreasing = FALSE)

    # get the bet two
    best_two_inds <- ordering[1:2]

    # take the average of the statistics
    best_two_average <- colMeans(rbind(hsn[best_two_inds[1],],
                                       hsn[best_two_inds[2],]))

    # now find the point that is GERGM_Object@convex_hull_proportion towards the
    # center:
    target_stats <- (best_two_average * GERGM_Object@convex_hull_proportion) +
      (simulated_center * (1 - GERGM_Object@convex_hull_proportion))

    # find the distance between the simulated center and the observed:
    obs_sim_center_distance <- calc_dists(observed_network_stats, simulated_center)

    # now see if the observed stats are closer to the center than 50% of the
    # simulated points
    closer_than_prop <- length(which(sim_center_distances > obs_sim_center_distance)) /
      length(sim_center_distances)

    cat("Distance between observed and mean simulated statistics:",
        obs_sim_center_distance,
        "\nProportion of simulated stats that are further away:",
        closer_than_prop,"\n")

    if (closer_than_prop > GERGM_Object@convex_hull_convergence_proportion) {
      converged <- TRUE
      cat("Converged!\n")
    } else {
      if (verbose) {
        theta.new <- optim(par = GERGM_Object@theta.par,
                           log.l,
                           alpha = GERGM_Object@weights,
                           hsnet = hsn,
                           ltheta = as.numeric(GERGM_Object@theta.par),
                           together = GERGM_Object@downweight_statistics_together,
                           possible.stats = possible.stats,
                           GERGM_Object = GERGM_Object,
                           override_statistics = target_stats, # select an easier target
                           method = GERGM_Object@optimization_method,
                           hessian = T,
                           control = list(fnscale = -1, trace = 6))
      } else {
        theta.new <- optim(par = GERGM_Object@theta.par,
                           log.l,
                           alpha = GERGM_Object@weights,
                           hsnet = hsn,
                           ltheta = as.numeric(GERGM_Object@theta.par),
                           together = GERGM_Object@downweight_statistics_together,
                           possible.stats = possible.stats,
                           GERGM_Object = GERGM_Object,
                           override_statistics = target_stats, # select an easier target
                           method = GERGM_Object@optimization_method,
                           hessian = T,
                           control = list(fnscale = -1, trace = 0))
      }

      # update the theta parameters.
      GERGM_Object@theta.par <- as.numeric(theta.new$par)

      cat("New thetas:",as.numeric(theta.new$par),"\n")
    }

  }




  # finally, return GREGM object:
  return(GERGM_Object)
}
