optimize_initialization <- function(GERGM_Object,
                                    verbose,
                                    seed2,
                                    possible.stats,
                                    theta,
                                    statistics){

  # now we do a grid search
  num_thetas <- length(theta$par)

  grid_steps  <- GERGM_Object@theta_grid_optimization_list$grid_steps
  step_size <- GERGM_Object@theta_grid_optimization_list$step_size
  total_steps <- 2 * grid_steps + 1
  cores <- GERGM_Object@theta_grid_optimization_list$cores

  multiplier <- GERGM_Object@theta_grid_optimization_list$iteration_fraction
  cur1 <- GERGM_Object@number_of_simulations
  GERGM_Object@number_of_simulations <- cur1 * multiplier
  cur2 <- GERGM_Object@burnin
  GERGM_Object@burnin <- cur2 * multiplier

  # generate the parameter sweep grid
  parameter_list <- vector(mode = "list", length = num_thetas)
  for (i in 1:num_thetas) {
    min_bnd <- theta$par[i] - grid_steps * step_size * abs(theta$par[i])
    max_bnd <- theta$par[i] + grid_steps * step_size * abs(theta$par[i])
    parameter_list[[i]] <- seq(from = min_bnd,
                               to = max_bnd,
                               length.out = total_steps)
  }

  parameter_grid <- data.frame(expand.grid(parameter_list))
  grid_size <- nrow(parameter_grid)

  vec <- 1:grid_size
  cat("Performing theta optimization grid search in parallel on",cores,
      "cores. Total grid size is",grid_size,
      "parameter combinations.This may take a while...\n")
  cl <- parallel::makeCluster(getOption("cl.cores", cores))

  results <- parallel::clusterApplyLB(cl = cl,
    x = vec,
    fun = theta_grid_search,
    parameter_grid = parameter_grid,
    GERGM_Object = GERGM_Object,
    seed2 = seed2,
    possible.stats = possible.stats,
    verbose = verbose,
    statistics = statistics)

  # stop the cluster when we are done
  parallel::stopCluster(cl)

  # trasnform into a vector
  differences <- as.numeric(unlist(results))

  # find the minimum difference
  min_diff <- which(differences == min(differences))

  optimal_thetas <- parameter_grid[min_diff,]

  cat("The minimum aggregate absolute difference was",min(differences),
      "for theta values",optimal_thetas,
      ". If this value is very large, a wider grid search may be necessary.\n")
  return(optimal_thetas)
} # end of function call
