
mple_distribution <- function(GERGM_Object,
                              verbose) {
  est <- GERGM_Object@theta.par
  ests <- NULL
  if (verbose) {
    ests <- optim(par = est,
                  pl_distribution,
                  GERGM_Object = GERGM_Object,
                  method = "BFGS",
                  hessian = TRUE,
                  control = list(fnscale = -1, trace = 6))
  } else {
    ests <- optim(par = est,
                  pl_distribution,
                  GERGM_Object = GERGM_Object,
                  method = "BFGS",
                  hessian = TRUE,
                  control = list(fnscale = -1, trace = 0))
  }
  return(ests)
}


# current version that works with all of our flexible new statistics
pl_distribution <- function(theta,
                            GERGM_Object){

  cat("Weighted MPLE Theta = ",theta,"\n")
  current_network <- GERGM_Object@network
  number_of_nodes <- nrow(current_network)
  triples <- GERGM_Object@statistic_auxiliary_data$triples
  pairs <- GERGM_Object@statistic_auxiliary_data$pairs

  dw <- as.numeric(GERGM_Object@downweight_statistics_together)

  integration_interval <- seq(from = 0,
                              to = 1,
                              length.out = GERGM_Object@integration_intervals)

  sad <- GERGM_Object@statistic_auxiliary_data

  # here we are only going to work with our actual thetas/statistics, not the full suite.
  selected_rows_matrix <- sad$specified_selected_rows_matrix
  rows_to_use <- sad$specified_rows_to_use
  base_statistics_to_save <- sad$specified_base_statistics_to_save
  base_statistic_alphas <- sad$specified_base_statistic_alphas
  num_non_base_statistics <- sum(GERGM_Object@non_base_statistic_indicator)

  objective <- mple_distribution_objective(
    number_of_nodes,
    GERGM_Object@stats_to_use - 1,
    current_network,
    triples - 1,
    pairs - 1,
    selected_rows_matrix - 1,
    rows_to_use - 1,
    base_statistics_to_save - 1,
    base_statistic_alphas,
    num_non_base_statistics,
    GERGM_Object@non_base_statistic_indicator,
    theta,
    GERGM_Object@weights,
    dw,
    integration_interval,
    GERGM_Object@parallel)

  # try some regularization with optional regularization weight
  objective <- objective - GERGM_Object@regularization_weight * sum(abs(theta)^2)
  cat("Calculation complete, objective is:",objective,"\n\n")
  return(objective)
}
