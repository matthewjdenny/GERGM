theta_grid_search <- function(x,
                              parameter_grid,
                              GERGM_Object,
                              seed2,
                              possible.stats,
                              verbose,
                              statistics){

  thetas <- as.numeric(parameter_grid[x,])
  GERGM_Object@theta.par <- thetas

  # now optimize the proposal variance if we are using Metropolis Hasings
  if (GERGM_Object@hyperparameter_optimization) {
    if (GERGM_Object@estimation_method == "Metropolis") {
      GERGM_Object@proposal_variance <- Optimize_Proposal_Variance(
        GERGM_Object = GERGM_Object,
        seed2 = seed2,
        possible.stats = possible.stats,
        verbose = verbose)
      cat("Proposal variance optimization complete! Proposal variance is:",
          GERGM_Object@proposal_variance,"\n",
          "--------- END HYPERPARAMETER OPTIMIZATION ---------",
          "\n\n")
    }
  }

  GERGM_Object <- Simulate_GERGM(GERGM_Object,
                                 seed1 = seed2,
                                 possible.stats = possible.stats,
                                 verbose = verbose)

  hsn <- GERGM_Object@MCMC_output$Statistics[,which(statistics == 1)]

  theta.new <- optim(par = thetas,
                     log.l,
                     alpha = GERGM_Object@reduced_weights,
                     hsnet = hsn,
                     ltheta = as.numeric(thetas),
                     together = GERGM_Object@downweight_statistics_together,
                     possible.stats = possible.stats,
                     GERGM_Object = GERGM_Object,
                     method = "BFGS",
                     hessian = T,
                     control = list(fnscale = -1, trace = 0))

  new_thetas <- theta.new$par
  # calculate absolute difference
  absolute_difference <- sum(abs(new_thetas - thetas))

  return(absolute_difference)
}
