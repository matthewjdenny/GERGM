#' @title A Function to simulate networks from a GERGM with given theta parameters.
#' @description Simulates networks from a GERGM for a given set of parameter values.
#'
#' @param formula A formula object that specifies which statistics the user would
#' like to include while simulating the network, and the network the user is
#' providing as the initial network. Currently, the following statistics can be
#' specified: c("out2stars", "in2stars", 	"ctriads", "mutual", "ttriads").
#' @param mutual The theta value provided for the reciprocity parameter, defaults
#' to 0. Only statistics for structural terms included in formula will be used.
#' @param ttriads The theta value provided for the transitive triads parameter,
#' defaults to 0. Only statistics for structural terms included in formula will
#' be used.
#' @param ctriads The theta value provided for the cyclic triads parameter,
#' defaults to 0. Only statistics for structural terms included in formula will
#' be used.
#' @param in2stars The theta value provided for the in 2-stars parameter,
#' defaults to 0. Only statistics for structural terms included in formula will
#' be used.
#' @param out2stars The theta value provided for the out 2-starts parameter,
#' defaults to 0. Only statistics for structural terms included in formula will
#' be used.
#' @param twostars The theta value provided for the undirected 2-starts parameter,
#' defaults to 0. Only statistics for structural terms included in formula will
#' be used.
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
#' @param omit_intercept_term Defualts to FALSE, can be set to TRUE if the user wishes to omit the model intercept term.
#' @param simulate_correlation_network Defaults to FALSE. Experimental.
#' @examples
#' set.seed(12345)
#' net <- matrix(runif(100),10,10)
#' diag(net) <- 0
#' colnames(net) <- rownames(net) <- letters[1:10]
#' formula <- net ~ ttriads + in2stars
#'
#' test <- simulate_networks(formula,
#'  ttriads = 0.6,
#'  in2stars = -0.8,
#'  network_is_directed = TRUE,
#'  simulation_method = "Metropolis",
#'  number_of_networks_to_simulate = 10000,
#'  thin = 1/10,
#'  proposal_variance = 0.5,
#'  downweight_statistics_together = TRUE,
#'  MCMC_burnin = 1000,
#'  seed = 456)
#' @return A list object containing simulated networks and parameters used to
#' specify the simulation. See the $MCMC_Output field for simulated networks.
#' @export
simulate_networks <- function(formula,
  mutual = 0,
  ttriads = 0,
  ctriads = 0,
  in2star2 = 0,
  out2star2 = 0,
  twostars = 0,
  simulation_method = c("Metropolis","Gibbs"),
  network_is_directed = c(TRUE, FALSE),
  number_of_networks_to_simulate = 500,
  thin = 1,
  proposal_variance = 0.1,
  downweight_statistics_together = TRUE,
  MCMC_burnin = 100,
  seed = 123,
  omit_intercept_term = FALSE,
  simulate_correlation_network = FALSE
){

  # This is the main function to estimate a GERGM model

  # hard coded possible stats
  possible_structural_terms <- c("out2stars", "in2stars", "ctriads", "mutual", "ttriads","edges")
  possible_structural_terms_undirected <- c("twostars", "ttriads")
  possible_covariate_terms <- c("absdiff", "nodecov", "nodefactor", "sender", "receiver", "intercept")
  possible_network_terms <- "netcov"
  # possible_transformations <- c("cauchy", "logcauchy", "gaussian", "lognormal")

  #check for an edges statistic
  form<- as.formula(formula)
  parsed <- deparse(form)
  if(length(parsed) > 1){
    parsed <- paste0(parsed, collapse = " ")
  }
  if(grepl("edges",parsed)){
    stop("You may not specify an edges statistic.")
  }

  # set logical values for whether we are using MPLE only, whether the network
  # is directed, and which estimation method we are using as well as the
  # transformation type
  network_is_directed <- network_is_directed[1] #default is TRUE
  simulation_method <- simulation_method[1] #default is Gibbs
  transformation_type <- "cauchy"
  normalization_type <- "division"

  # check terms for undirected network
  if(!network_is_directed){
    formula <- parse_undirected_structural_terms(
      formula,
      possible_structural_terms,
      possible_structural_terms_undirected)
  }

  # automatically add an intercept term unless omit_intercept_term is TRUE
  if(!omit_intercept_term){
    formula <- add_intercept_term(formula)
  }

  # if we are using a correlation network, then the network must be undirected.
  if(simulate_correlation_network){
    network_is_directed <- FALSE
  }

  # convert logical to numeric indicator
  if(downweight_statistics_together){
    downweight_statistics_together <- 1
  }else{
    downweight_statistics_together <- 0
  }

  #make sure proposal variance is greater than zero
  #make sure proposal variance is greater than zero
  if(proposal_variance <= 0.001){
    proposal_variance <- 0.001
    cat("You supplied a proposal variance that was less than or equal to zero.
        It has been reset to 0.001, considder respecifying...\n")
  }

  formula <- as.formula(formula)

  #0. Prepare the data
  Transformed_Data <- Prepare_Network_and_Covariates(
    formula,
    possible_structural_terms,
    possible_covariate_terms,
    possible_network_terms,
    covariate_data = NULL,
    normalization_type = normalization_type,
    is_correlation_network = simulate_correlation_network,
    is_directed = network_is_directed)

  # create theta coefficients
  theta_coeficients = NULL
  if(out2stars != 0){
    theta_coeficients <- c(theta_coeficients, out2stars)
  }
  if(in2stars != 0 | twostars != 0){
    theta_coeficients <- c(theta_coeficients, in2stars)
  }
  if(ctriads != 0){
    theta_coeficients <- c(theta_coeficients, ctriads)
  }
  if(mutual != 0){
    theta_coeficients <- c(theta_coeficients, mutual)
  }
  if(ttriads != 0){
    theta_coeficients <- c(theta_coeficients, ttriads)
  }

  #1. Create GERGM object from network

  GERGM_Object <- Create_GERGM_Object_From_Formula(
    formula,
    theta.coef = theta_coeficients,
    possible_structural_terms,
    possible_covariate_terms,
    possible_network_terms,
    raw_network = Transformed_Data$network,
    together = 1,
    transform.data = NULL,
    lambda.coef = NULL,
    transformation_type = transformation_type,
    is_correlation_network = simulate_correlation_network,
    is_directed = network_is_directed
  )

  print(GERGM_Object@theta.coef)
  print(GERGM_Object@theta.par)

  GERGM_Object@theta_estimation_converged <- TRUE
  GERGM_Object@lambda_estimation_converged <- TRUE
  GERGM_Object@observed_network  <- GERGM_Object@network
  GERGM_Object@observed_bounded_network <- GERGM_Object@bounded.network
  GERGM_Object@simulation_only <- TRUE
  GERGM_Object@theta.par <- theta_coeficients

  if(network_is_directed){
    GERGM_Object@undirected_network <- FALSE
  }else{
    GERGM_Object@undirected_network <- TRUE
  }


  # if we are using a correlation network then set field to TRUE.
  if(simulate_correlation_network){
    GERGM_Object@is_correlation_network <- TRUE
  }else{
    GERGM_Object@is_correlation_network <- FALSE
  }


  #now simulate from last update of theta parameters
  GERGM_Object <- Simulate_GERGM(GERGM_Object,
                                 nsim = number_of_networks_to_simulate,
                                 method = simulation_method,
                                 MCMC.burnin = MCMC_burnin,
                                 thin = thin,
                                 shape.parameter = proposal_variance,
                                 together = downweight_statistics_together,
                                 seed1 = seed,
                                 possible.stats = possible_structural_terms)

  #which(GERGM_Object@stats_to_use == 1)
  num.nodes <- GERGM_Object@num_nodes
  triples <- t(combn(1:num.nodes, 3))
  pairs <- t(combn(1:num.nodes, 2))
  # initialize the network with the observed network
  init.statistics <- h2(GERGM_Object@bounded.network,
                        triples = triples,
                        statistics = rep(1, length(possible_structural_terms)),
                        alphas = GERGM_Object@weights,
                        together = downweight_statistics_together)

  hsn.tot <- GERGM_Object@MCMC_output$Statistics
  #calculate t.test p-values for calculating the difference in the means of


  stats.data <- data.frame(Observed = init.statistics,
                           Simulated = colMeans(hsn.tot))
  rownames(stats.data) <- possible_structural_terms
  cat("Statistics of initial network and networks simulated from given theta
      parameter estimates:\n")

  print(stats.data)

  # change back column names if we are dealing with an undirected network
  if(!network_is_directed){
    change <- which(colnames(GERGM_Object@theta.coef) == "in2star")
    if(length(change) > 0){
      colnames(GERGM_Object@theta.coef)[change] <- "twostars"
    }
  }

  # make GOF plot

  GOF(GERGM_Object)
  Sys.sleep(2)
  Trace_Plot(GERGM_Object)

  return_list <- list(formula = GERGM_Object@formula,
                      theta_values = GERGM_Object@theta.coef[,1],
                      alpha_values = GERGM_Object@weights,
                      MCMC_Output = GERGM_Object@MCMC_output)
  return(return_list)
}
