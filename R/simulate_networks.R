#' @title A Function to simulate networks from a GERGM with given theta parameters.
#' @description Simulates networks from a GERGM for a given set of parameter values.
#'
#' @param formula A formula object that specifies which statistics the user would
#' like to include while simulating the network, and the network the user is
#' providing as the initial network. Currently, the following statistics can be
#' specified: c("edges", "out2stars", "in2stars", 	"ctriads", "mutual",
#' "ttriads").
#' @param thetas A vector of theta parameters given in the same order as the
#' formula terms, which the user would like to use to parameterize the model.
#' @param network_is_directed Logical specifying whether or not the observed
#' network is directed. Default is TRUE.
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
#' @param downweight_statistics_together Logical specifying whether or not the
#' weights should be applied inside or outside the sum. Default is TRUE and user
#' should not select FALSE under normal circumstances.
#' @param MCMC_burnin Number of samples from the MCMC simulation procedure that
#' will be discarded before drawing the samples used for estimation. Default is
#' 100.
#' @param seed Seed used for reproducibility. Default is 123.
#' @param GERGM_Object Optional argument allowing the user to supply a GERGM
#' object output by the gergm() estimation function in order to simualt further
#' networks. Defaults to NULL. If a GERGM object is provided, any user specified
#' parameter values will be ignored and the final parameter estimates from the
#' gergm() function will be used instead. When using this option, the folowing
#' terms must still be specified: number_of_networks_to_simulate, thin, and
#' MCMC_burnin. proposal_variance may also be specified, or if set equal to NULL,
#' then the proposal variance from parameter estimation will be used instead (
#' this option is likely preferred in most situations).
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
#' @param covariate_data A data frame containing node level covariates the user
#' wished to transform into sender or reciever effects. It must have row names
#' that match every entry in colnames(raw_network), should have descriptive
#' column names.  If left NULL, then no sender or reciever effects will be
#' added.
#' @param lambdas A vector of lambda parameters given in the same order as the
#' formula terms, which the user would like to use to parameterize the model.
#' Covariate effects should be specified after endogenous effects.
#' @param ... Optional arguments, currently unsupported.
#' @examples
#' set.seed(12345)
#' net <- matrix(runif(100),10,10)
#' diag(net) <- 0
#' colnames(net) <- rownames(net) <- letters[1:10]
#' formula <- net ~ edges + ttriads + in2stars
#'
#' test <- simulate_networks(formula,
#'  thetas = c(0.6,-0.8),
#'  lambdas = 0.2,
#'  number_of_networks_to_simulate = 100,
#'  thin = 1/10,
#'  proposal_variance = 0.5,
#'  MCMC_burnin = 1000,
#'  seed = 456)
#'
#' # preferred method for specifying a null model
#' formula <- net ~ edges(method = "endogenous")
#' test <- simulate_networks(
#'  formula,
#'  thetas = 0,
#'  number_of_networks_to_simulate = 1000,
#'  thin = 1/10,
#'  proposal_variance = 0.5,
#'  MCMC_burnin = 1000,
#'  seed = 456)
#' @return A list object containing simulated networks and parameters used to
#' specify the simulation. See the $MCMC_Output field for simulated networks. If
#' GERGM_Object is provded, then a GERGM object will be returned instead.
#' @export
simulate_networks <- function(formula,
  thetas,
  simulation_method = c("Metropolis","Gibbs"),
  network_is_directed = TRUE,
  number_of_networks_to_simulate = 500,
  thin = 1,
  proposal_variance = 0.1,
  downweight_statistics_together = TRUE,
  MCMC_burnin = 100,
  seed = 123,
  GERGM_Object = NULL,
  return_constrained_networks = FALSE,
  optimize_proposal_variance = FALSE,
  target_accept_rate = 0.25,
  use_stochastic_MH = FALSE,
  stochastic_MH_proportion = 1,
  covariate_data = NULL,
  lambdas = NULL,
  ...
){

  GERGM_Object_Provided <- TRUE
  if (is.null(GERGM_Object)) {
    GERGM_Object_Provided <- FALSE
  }

  # hard coded possible stats
  possible_structural_terms <- c("out2stars", "in2stars", "ctriads", "mutual", "ttriads","edges")
  possible_structural_terms_undirected <- c("twostars",
                                            "ttriads",
                                            "edges")
  if (network_is_directed) {
    possible_structural_term_indices <- 1:6
  } else {
    possible_structural_term_indices <- c(2,5,6)
  }
  possible_covariate_terms <- c("absdiff", "nodecov", "nodematch", "sender", "receiver", "intercept", "nodemix")
  possible_network_terms <- "netcov"
  # possible_transformations <- c("cauchy", "logcauchy", "gaussian", "lognormal")

  if (is.null(GERGM_Object)) {
    # pass in experimental correlation network feature through elipsis
    weighted_MPLE <- FALSE

    # This is the main function to estimate a GERGM model

    # check to see if we are using a regression intercept term, and if we are,
    # then add the "intercept" field to the formula.
    check_for_edges <- Parse_Formula_Object(formula,
                                            possible_structural_terms,
                                            possible_covariate_terms,
                                            possible_network_terms,
                                            raw_network = NULL,
                                            theta = NULL,
                                            terms_to_parse = "structural",
                                            covariate_data = covariate_data)

    # automatically add an intercept term if necessary
    if (check_for_edges$include_intercept) {
      if (check_for_edges$regression_intercept) {
        formula <- add_intercept_term(formula)
      }
    }

    # set logical values for whether we are using MPLE only, whether the network
    # is directed, and which estimation method we are using as well as the
    # transformation type
    network_is_directed <- network_is_directed[1] #default is TRUE
    simulation_method <- simulation_method[1] #default is Gibbs
    transformation_type <- "cauchy"
    normalization_type <- "division"

    # check terms for undirected network
    if (!network_is_directed) {
      formula <- parse_undirected_structural_terms(
        formula,
        possible_structural_terms,
        possible_structural_terms_undirected)
    }

    #make sure proposal variance is greater than zero
    if (proposal_variance <= 0) {
      stop("You supplied a proposal variance that was less than or equal to zero.")
    }

    formula <- as.formula(formula)

    #0. Prepare the data
    Transformed_Data <- Prepare_Network_and_Covariates(
      formula,
      possible_structural_terms,
      possible_covariate_terms,
      possible_network_terms,
      covariate_data = covariate_data,
      normalization_type = normalization_type,
      is_directed = network_is_directed)

    data_transformation <- NULL
    if (!is.null(Transformed_Data$transformed_covariates)) {
      data_transformation <- Transformed_Data$transformed_covariates
    }

    #1. Create GERGM object from network
    GERGM_Object <- Create_GERGM_Object_From_Formula(
      formula,
      theta.coef = thetas,
      possible_structural_terms,
      possible_covariate_terms,
      possible_network_terms,
      raw_network = Transformed_Data$network,
      together = 1,
      transform.data = data_transformation,
      lambda.coef = lambdas,
      transformation_type = transformation_type,
      is_directed = network_is_directed,
      possible_structural_terms_undirected = possible_structural_terms_undirected)

    if (!is.null(data_transformation)) {
      GERGM_Object@data_transformation <- data_transformation
    }

    GERGM_Object@transformation_type <- transformation_type
    GERGM_Object@theta_estimation_converged <- TRUE
    GERGM_Object@lambda_estimation_converged <- TRUE
    GERGM_Object@observed_network  <- GERGM_Object@network
    GERGM_Object@observed_bounded_network <- GERGM_Object@bounded.network
    GERGM_Object@simulation_only <- TRUE
    GERGM_Object@theta.par <- thetas
    GERGM_Object@directed_network <- network_is_directed
    GERGM_Object@proposal_variance <- proposal_variance
    GERGM_Object@estimation_method <- simulation_method
    GERGM_Object@target_accept_rate <- target_accept_rate
    GERGM_Object@number_of_simulations <- number_of_networks_to_simulate
    GERGM_Object@thin <- thin
    GERGM_Object@burnin <- MCMC_burnin
    GERGM_Object@downweight_statistics_together <- downweight_statistics_together
    GERGM_Object@weighted_MPLE <- weighted_MPLE
    GERGM_Object@hyperparameter_optimization <- optimize_proposal_variance
    GERGM_Object@fine_grained_pv_optimization <- TRUE
    GERGM_Object@parallel <- FALSE
    GERGM_Object@parallel_statistic_calculation <- FALSE
    GERGM_Object@cores <- 1
    GERGM_Object@use_stochastic_MH <- use_stochastic_MH
    GERGM_Object@stochastic_MH_proportion <- stochastic_MH_proportion
    GERGM_Object@possible_endogenous_statistic_indices <- possible_structural_term_indices

    # prepare auxiliary data
    GERGM_Object@statistic_auxiliary_data <- prepare_statistic_auxiliary_data(
      GERGM_Object)

    # record statistics on observed and bounded network
    h.statistics1 <- calculate_h_statistics(
      GERGM_Object,
      GERGM_Object@statistic_auxiliary_data,
      all_weights_are_one = FALSE,
      calculate_all_statistics = TRUE,
      use_constrained_network = FALSE)
    h.statistics2 <- calculate_h_statistics(
      GERGM_Object,
      GERGM_Object@statistic_auxiliary_data,
      all_weights_are_one = FALSE,
      calculate_all_statistics = TRUE,
      use_constrained_network = TRUE)
    statistic_values <- rbind(h.statistics1, h.statistics2)
    colnames(statistic_values) <- GERGM_Object@full_theta_names
    rownames(statistic_values) <- c("network", "bounded_network")
    GERGM_Object@stats <- statistic_values

    # allow the user to specify mu and phi
     if (check_for_edges$regression_intercept) {
      beta <- lambdas[1:(length(lambdas) - 1)]
      sig <- 0.01 + exp(lambdas[length(lambdas)])
      BZ <- 0
      for (j in 1:(dim(GERGM_Object@data_transformation)[3])) {
        BZ <- BZ + beta[j] * GERGM_Object@data_transformation[, , j]
      }

      #store so we can transform back
      GERGM_Object@BZ <- BZ
      GERGM_Object@BZstdev <- sig
    }

  } else {
    # a GERGM_Object was provided
    if (!is.null(proposal_variance)) {
      GERGM_Object@proposal_variance <- proposal_variance
    }
    GERGM_Object@number_of_simulations <- number_of_networks_to_simulate
    GERGM_Object@thin <- thin
    GERGM_Object@burnin <- MCMC_burnin
    network_is_directed <- GERGM_Object@directed_network
  }


  if (GERGM_Object@hyperparameter_optimization){
    if (GERGM_Object@estimation_method == "Metropolis") {
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

  #now simulate from last update of theta parameters
  GERGM_Object <- Simulate_GERGM(GERGM_Object,
                                 seed1 = seed,
                                 possible.stats = possible_structural_terms)

  # initialize the network with the observed network
  init.statistics <- calculate_h_statistics(
    GERGM_Object,
    GERGM_Object@statistic_auxiliary_data,
    all_weights_are_one = FALSE,
    calculate_all_statistics = TRUE,
    use_constrained_network = TRUE)

  hsn.tot <- GERGM_Object@MCMC_output$Statistics

  # save these statistics so we can make GOF plots in the future, otherwise
  # they would be the transformed statistics which would produce poor GOF plots.
  GERGM_Object@simulated_statistics_for_GOF <- hsn.tot


  stats.data <- data.frame(Observed = init.statistics,
                           Simulated = colMeans(hsn.tot))
  rownames(stats.data) <- GERGM_Object@full_theta_names
  cat("Statistics of initial network and networks simulated from given theta
      parameter estimates:\n")

  print(stats.data)

  # fix constrained values
  GERGM_Object@stats[2,] <- init.statistics

  # change back column names if we are dealing with an undirected network
  if (!network_is_directed) {
    change <- which(colnames(GERGM_Object@theta.coef) == "in2star")
    if (length(change) > 0) {
      colnames(GERGM_Object@theta.coef)[change] <- "twostars"
    }
  }

  GERGM_Object@simulated_bounded_networks_for_GOF <- GERGM_Object@MCMC_output$Networks

  # precalculate intensities and degree distributions so we can return them.
  GERGM_Object <- calculate_additional_GOF_statistics(GERGM_Object)


  # make GOF plot
  try({
    GOF(GERGM_Object)
    Sys.sleep(2)
    Trace_Plot(GERGM_Object)
  })

  if (!return_constrained_networks) {
    cat("Transforming networks simulated via MCMC as part of the fit diagnostics back on to the scale of observed network. You can access these networks through the '@MCMC_output$Networks' field returned by this function...\n")
    GERGM_Object@simulated_bounded_networks_for_GOF <- GERGM_Object@MCMC_output$Networks
    GERGM_Object <- Convert_Simulated_Networks_To_Observed_Scale(GERGM_Object)
  } else {
    cat("Returning constrained [0,1] simulated networks...\n")
  }

  if (!GERGM_Object_Provided) {
    return_list <- list(formula = GERGM_Object@formula,
                        theta_values = GERGM_Object@theta.coef[,1],
                        alpha_values = GERGM_Object@weights,
                        MCMC_Output = GERGM_Object@MCMC_output,
                        bounded_networks = GERGM_Object@simulated_bounded_networks_for_GOF)
    return(return_list)
  } else {
    return(GERGM_Object)
  }

}
