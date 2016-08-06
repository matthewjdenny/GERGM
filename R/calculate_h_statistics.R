# statistic_auxiliary_data <- prepare_statistic_auxiliary_data(GERGM_Object)
# calculate_h_statistics(GERGM_Object, statistic_auxiliary_data,calculate_all_statistics = TRUE)
calculate_h_statistics <- function(GERGM_Object,
                                   statistic_auxiliary_data,
                                   all_weights_are_one = FALSE,
                                   calculate_all_statistics = FALSE,
                                   use_constrained_network = TRUE,
                                   network = NULL) {

  # extract parameters to be used by the c++ function
  num_non_base_statistics <- sum(GERGM_Object@non_base_statistic_indicator)
  num_nodes <- GERGM_Object@num_nodes
  triples <- t(combn(1:num_nodes, 3))
  pairs <- t(combn(1:num_nodes, 2))

  alphas <- GERGM_Object@weights
  if (all_weights_are_one) {
    alphas <- rep(1, alphas)
  }
  together <- as.numeric(GERGM_Object@downweight_statistics_together)

  # if no network was provided, get it.
  if (is.null(network)) {
    if (use_constrained_network) {
      network <- GERGM_Object@bounded.network
    } else {
      network <- GERGM_Object@network
    }
  }

  if (calculate_all_statistics) {
    selected_rows_matrix <- statistic_auxiliary_data$full_selected_rows_matrix
    rows_to_use <- statistic_auxiliary_data$specified_rows_to_use
    base_statistics_to_save <- statistic_auxiliary_data$full_base_statistics_to_save
    base_statistic_alphas <- statistic_auxiliary_data$full_base_statistic_alphas
  } else {
    selected_rows_matrix <- statistic_auxiliary_data$specified_selected_rows_matrix
    rows_to_use <- statistic_auxiliary_data$specified_rows_to_use
    base_statistics_to_save <- statistic_auxiliary_data$specified_base_statistics_to_save
    base_statistic_alphas <- statistic_auxiliary_data$specified_base_statistic_alphas
  }
  h_stats <- h_statistics(
    GERGM_Object@stats_to_use - 1,
    network,
    triples - 1,
    pairs - 1,
    alphas,
    together,
    selected_rows_matrix - 1,
    rows_to_use - 1,
    base_statistics_to_save - 1,
    base_statistic_alphas,
    num_non_base_statistics,
    GERGM_Object@non_base_statistic_indicator)

  # now if we are not returning all of the statistics, then reorder the
  # statistics we are returning so they match up with the statistics_to_use
  # ordering
  if (!calculate_all_statistics) {
    reordering <- c(which(GERGM_Object@non_base_statistic_indicator == 0),
                    which(GERGM_Object@non_base_statistic_indicator == 1))
    h_stats <- h_stats[reordering]
  }
  return(as.numeric(h_stats))
}
