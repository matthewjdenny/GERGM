#' @title A Function to predict edge weights from a GERGM fit object.
#' @description Performs edgewise predictions from a GERGM model fit.
#'
#' @param GERGM_Object A GERGM object output by the gergm() estimation function.
#' The folowing terms must still be specified: number_of_networks_to_simulate,
#' thin, and MCMC_burnin. proposal_variance may also be specified, or if set
#' equal to NULL, then the proposal variance from parameter estimation will be
#' instead (this option is likely preferred in most situations).
#' @param simulation_method Default is "Metropolis" which allows for exponential
#' downweighting, can also be "Gibbs".
#' @param number_of_networks_to_simulate Number of simulations generated for
#' estimation via MCMC. Default is 500.
#' @param thin The proportion of samples that are kept from each simulation. For
#' example, thin = 1/200 will keep every 200th network in the overall simulated
#' sample. Default is 1.
#' @param proposal_variance The variance specified for the Metropolis Hastings
#' simulation method. This parameter is inversely proportional to the average
#' acceptance rate of the M-H sampler and should be adjusted so that the average
#' acceptance rate is approximately 0.25. Default is 0.1.
#' @param MCMC_burnin Number of samples from the MCMC simulation procedure that
#' will be discarded before drawing the samples used for estimation. Default is
#' 100.
#' @param seed Seed used for reproducibility. Default is 123.
#' @param return_constrained_networks Logical argument indicating whether
#' simulated networks should be transformed back to observed scale or whether
#' constrained [0,1] networks should be returned. Defaults to FALSE, in which
#' case networks are returned on observed scale.
#' @param optimize_proposal_variance Logical indicating whether proposal
#' variance should be optimized if using Metropolis Hastings for simulation.
#' Defaults to FALSE.
#' @param target_accept_rate Defaults to 0.25, can be used to optimize
#' Metropolis Hastings simulations.
#' @param use_stochastic_MH A logical indicating whether a stochastic approximation
#' to the h statistics should be used under Metropolis Hastings in-between
#' thinned samples. This may dramatically speed up estimation. Defualts to FALSE.
#' HIGHLY EXPERIMENTAL!
#' @param stochastic_MH_proportion Percentage of dyads/triads to use for
#' approximation, defaults to 0.25
#' @return A list object containing simulated networks.
#' @export
conditional_edge_prediction <- function(
  GERGM_Object,
  simulation_method = c("Metropolis","Gibbs"),
  number_of_networks_to_simulate = 500,
  thin = 1,
  proposal_variance = 0.1,
  MCMC_burnin = 100,
  seed = 123,
  return_constrained_networks = FALSE,
  optimize_proposal_variance = FALSE,
  target_accept_rate = 0.25,
  use_stochastic_MH = FALSE,
  stochastic_MH_proportion = 1
){



  # a GERGM_Object was provided
  if (!is.null(proposal_variance)) {
    GERGM_Object@proposal_variance <- proposal_variance
  }
  GERGM_Object@number_of_simulations <- number_of_networks_to_simulate
  GERGM_Object@thin <- thin
  GERGM_Object@burnin <- MCMC_burnin
  network_is_directed <- GERGM_Object@directed_network

  # hard coded possible stats
  possible_structural_terms <- c("out2stars", "in2stars", "ctriads", "mutual", "ttriads","edges")

  if (network_is_directed) {
    possible_structural_term_indices <- 1:6
  } else {
    possible_structural_term_indices <- c(2,5,6)
  }

  if (optimize_proposal_variance){
    if (simulation_method == "Metropolis") {
      GERGM_Object@proposal_variance <- Optimize_Proposal_Variance(
        GERGM_Object = GERGM_Object,
        seed2 = seed,
        possible.stats = possible_structural_terms,
        verbose = TRUE,
        max_updates = 50,
        fine_grained_optimization = TRUE,
        iteration_fraction = 1)
      cat("Proposal variance optimization complete! Proposal variance is:",
          GERGM_Object@proposal_variance,"\n",
          "--------- END HYPERPARAMETER OPTIMIZATION ---------",
          "\n\n")
    }
  }


  Edge_Predictions <- edge_predictions_bounded_scale <- NULL
  num_combinations <- GERGM_Object@num_nodes * (GERGM_Object@num_nodes -1)
  counter <- 1
  # get samples of each edge value and store them in an array.
  for (i in 1:GERGM_Object@num_nodes) {
    for (j in 1:GERGM_Object@num_nodes) {
      if (i != j) {
        cat("Predicting edge:",counter, "of",num_combinations,"...\n")
        counter <- counter + 1
        #now simulate from last update of theta parameters
        temp <- Simulate_GERGM(GERGM_Object,
                                       seed1 = seed,
                                       possible.stats = possible_structural_terms,
                                       predict_conditional_edges = TRUE,
                                       i = i,
                                       j = j)

        simulated_scale <- temp@MCMC_output$Networks[i,j,]

        # covert back to observed scale (so we incorporate in covariate effects)
        temp <- Convert_Simulated_Networks_To_Observed_Scale(temp)

        #now get and save the edge sample
        edge_values <- temp@MCMC_output$Networks[i,j,]

        # initialize this way to make sure we get dimensions right
        if (is.null(Edge_Predictions)) {
          
            Edge_Predictions <- array(0,dim = c(GERGM_Object@num_nodes,
                                                GERGM_Object@num_nodes,
                                                length(edge_values)))
            edge_predictions_bounded_scale <- array(0,dim = c(GERGM_Object@num_nodes,
                                                              GERGM_Object@num_nodes,
                                                              length(edge_values)))
          
        }

        # fill in the spot in the array
        Edge_Predictions[i,j,] <- edge_values
        edge_predictions_bounded_scale[i,j,] <- simulated_scale
      }
    }
  }

  sim_means <- apply(edge_predictions_bounded_scale,c(1,2),mean)
  sim_sds <- apply(edge_predictions_bounded_scale,c(1,2),sd)
  means <- apply(Edge_Predictions,c(1,2),mean)
  sds <- apply(Edge_Predictions,c(1,2),sd)

  # calculate the precentage of edges in CI
  in_ci <- rep(0,num_combinations)
  average_ci_coverage <- rep(0,num_combinations)
  average_sim_ci_coverage <- rep(0,num_combinations)
  counter <- 1
  for (i in 1:GERGM_Object@num_nodes) {
    for (j in 1:GERGM_Object@num_nodes) {
      if (i != j) {
        cur <- GERGM_Object@network[i,j]
        lower <- means[i,j] - 1.96 * sds[i,j]
        upper <- means[i,j] + 1.96 * sds[i,j]
        sim_lower <- sim_means[i,j] - 1.96 * sim_sds[i,j]
        sim_upper <- sim_means[i,j] + 1.96 * sim_sds[i,j]
        if(lower < cur  & cur < upper){
          in_ci[counter] <- 1
        }
        average_ci_coverage[counter] <- upper - lower
        average_sim_ci_coverage[counter] <- sim_upper - sim_lower
        counter <- counter + 1
      }
    }
  }


  cat("\n\n\n===========================================\n\nThe",
    "proportion of observed edges in predicted",
    "edge 95% confidence intervals is:", mean(in_ci),"\n\n")

  
    cat("Average confidence interval coverage is ",
        round(100*mean(average_sim_ci_coverage),2),
        "% of the possible simulated network space... \n",sep = "")
  


  # now lets do some comparison work.
  # GERGM_Object@bounded.network

  return(list(edge_predictions = Edge_Predictions,
              observed_network = GERGM_Object@network,
              predicted_edge_values = means,
              predicted_edge_sds = sds,
              proportion_in_ci = mean(in_ci),
              average_ci_coverage = 100*mean(average_ci_coverage),
              edge_predictions_bounded_scale = edge_predictions_bounded_scale))
}