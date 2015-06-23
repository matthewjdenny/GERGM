#' A Function to estimate a GERGM.
#'
#' @param formula A formula object that specifies the relationship between statistics and the observed network. Currently, the following statistics can be specified: c("out2star", "in2star", 	"ctriads", "recip", "ttriads", "edgeweight").
#' @param network_is_directed Logical specifying whether or not the observed network is directed. Default is TRUE.
#' @param use_MPLE_only Logical specifying whether or not only the maximum pseudo likelihood estimates should be obtained. In this case, no simulations will be performed. Default is FALSE.
#' @param data_transformation An n x n x m array where each of m layers contains a covariate that models the transform of the unbounded weighted network to a network whose edges are all on the unit interval. Default is NULL.
#' @param estimation_method Simulation method for MCMC estimation. Default is "Gibbs" which will generally be faster with well behaved networks but will not allow for exponential downweighting.
#' @param maximum_number_of_lambda_updates Maximum number of iterations of outer MCMC loop which alternately estimates transform parameters and ERGM parameters. In the case that data_transformation = NULL, this argument does not matter. Default is 10.
#' @param maximum_number_of_theta_updates Maximum number of iterations within the MCMC inner loop which estimates the ERGM parameters. Default is 100.
#' @param number_of_networks_to_simulate Number of simulations generated for estimation via MCMC. Default is 500.
#' @param thin The proportion of samples that are kept from each simulation. For example, thin = 1/200 will keep every 200th network in the overall simulated sample. Default is 1.
#' @param proposal_variance The variance specified for the Metropolis Hastings simulation method. This parameter is inversely proportional to the average acceptance rate of the M-H sampler and should be adjusted so that the average acceptance rate is approximately 0.25. 		Default is 0.1.
#' @param exponential_weights A vector of weights specifying the down weighting (via exponentiation) of each possible statistic. This vector must be the same length as the number of statistics used in object. Values are between 0 and 1. Default is NULL specifying a 1 for each statistic.
#' @param downweight_statistics_together Logical specifying whether or not the weights should be applied inside or outside the sum. Default is TRUE and user should not select FALSE under normal circumstances.
#' @param MCMC_burnin Number of samples from the MCMC simulation procedure that will be discarded before drawing the samples used for estimation. Default is 100.
#' @param seed Seed used for reproducibility. Default is 123.
#' @param convergence_tolerance Threshold designated for stopping criterion. If the difference of parameter estimates from one iteration to the next all have a p-value (under a paired t-test) greater than this value, the parameter estimates are declared to have converged. Default is 0.01.
#' @param MPLE_gain_factor Multiplicative constant between 0 and 1 that controls how far away the initial theta estimates will be from the standard MPLEs via a one step Fisher update. In the case of strongly dependent data, it is suggested to use a value of 0.10. Default is 0.
#' @param acceptable_fit_p_value_threshold A p-value threshold for how closely statistics of observed network conform to statistics of networks simulated from GERGM parameterized by converged final parameter estimates. Default value is 0.05.
#' @param force_x_theta_updates Defaults to 1 where theta estimation is not allowed to converge until thetas have updated for x iterations . Useful when model is not degenerate but simulated statistics do not match observed network well when algorithm stops after first x updates.
#' @return A gergm object containing parameter estimates.
#' @export
gergm <- function(formula,
                  network_is_directed = c(TRUE, FALSE),
                  use_MPLE_only = c(FALSE, TRUE),
                  data_transformation = NULL,
                  estimation_method = c("Gibbs", "Metropolis"),
                  maximum_number_of_lambda_updates = 10,
                  maximum_number_of_theta_updates = 100,
                  number_of_networks_to_simulate = 500,
                  thin = 1,
                  proposal_variance = 0.1,
                  exponential_weights = NULL,
                  downweight_statistics_together = TRUE,
                  MCMC_burnin = 100,
                  seed = 123,
                  convergence_tolerance = 0.01,
                  MPLE_gain_factor = 0,
                  acceptable_fit_p_value_threshold = 0.05,
                  force_x_theta_updates = 1){

  #' This is the main function to estimate a GERGM model

  #' hard coded possible stats
  possible.stats <- c("out2star", "in2star", "ctriads", "recip", "ttriads",
                      "edgeweight")

  #' set logical values for whether we are using MPLE only, whether the network
  #' is directed, and which estimation method we are using
  use_MPLE_only <- use_MPLE_only[1] #default is FALSE
  network_is_directed <- network_is_directed[1] #default is TRUE
  estimation_method <- estimation_method[1] #default is Gibbs

  #' convert logical to numeric indicator
  if(downweight_statistics_together){
    downweight_statistics_together <- 1
  }else{
    downweight_statistics_together <- 0
  }

  formula <- as.formula(formula)

  #1. Create GERGM object from network

  GERGM_Object <- Create_GERGM_Object_From_Formula(formula,
                                                   theta.coef = NULL,
                                                   possible.stats,
                                                   together = 1,
                                                   weights = exponential_weights,
                                                   transform.data = data_transformation,
                                                   lambda.coef = NULL)

  GERGM_Object@theta_estimation_converged <- FALSE
  GERGM_Object@lambda_estimation_converged <- FALSE
  GERGM_Object@observed_network  <- GERGM_Object@network
  GERGM_Object@observed_bounded_network <- GERGM_Object@bounded.network
  if(!is.null(data_transformation)){
    GERGM_Object@data_transformation <- data_transformation
  }


  #2. Estimate GERGM
  GERGM_Object <- Estimate_GERGM(formula,
                                 directed = network_is_directed,
                                 MPLE.only = use_MPLE_only,
                                 transform.data = data_transformation,
                                 method = estimation_method,
                                 max.num.iterations = maximum_number_of_lambda_updates,
                                 mc.num.iterations = maximum_number_of_theta_updates,
                                 nsim = number_of_networks_to_simulate,
                                 thin = thin,
                                 shape.parameter = proposal_variance,
                                 exponential_weights = exponential_weights,
                                 together = downweight_statistics_together,
                                 MCMC.burnin = MCMC_burnin,
                                 seed = seed,
                                 tolerance = convergence_tolerance,
                                 gain.factor = MPLE_gain_factor,
                                 possible.stats = possible.stats,
                                 GERGM_Object = GERGM_Object,
                                 force_x_theta_updates = force_x_theta_updates)

  #3. Perform degeneracy diagnostics and create GOF plots
  if(!GERGM_Object@theta_estimation_converged){
    warning("Estimation proceedure did not detect convergence in Theta estimates. Estimation halted when maximum number of updates was reached. Be careful to assure good model fit or select a more relaxed convergence criterion.")
  }
  if(!GERGM_Object@lambda_estimation_converged){
    warning("Estimation proceedure did not detect convergence in Lambda estimates. Estimation halted when maximum number of updates was reached. Be careful to assure good model fit or select a more relaxed convergence criterion.")
  }

  #now simulate from last update of theta parameters and
  GERGM_Object <- Simulate_GERGM(GERGM_Object,
                                 nsim = number_of_networks_to_simulate,
                                 method = estimation_method,
                                 MCMC.burnin = MCMC_burnin,
                                 thin = thin,
                                 shape.parameter = proposal_variance,
                                 together = downweight_statistics_together,
                                 seed1 = seed,
                                 possible.stats = possible.stats)

  #which(GERGM_Object@stats_to_use == 1)

  num.nodes <- GERGM_Object@num_nodes
  triples <- t(combn(1:num.nodes, 3))
  pairs <- t(combn(1:num.nodes, 2))
  # initialize the network with the observed network
  init.statistics <- h2(GERGM_Object@bounded.network,
                        triples = triples,
                        statistics = rep(1, length(possible.stats)),
                        alphas = GERGM_Object@weights,
                        together = downweight_statistics_together)

  hsn.tot <- GERGM_Object@MCMC_output$Statistics
  #calculate t.test p-values for calculating the difference in the means of
  # the newly simulated data with the original network
  statistic_test_p_values <- rep(NA,length(possible.stats))
  for(i in 1:length(possible.stats)){
    statistic_test_p_values[i] <- t.test(hsn.tot[, i],
                                      mu = init.statistics[i])$p.value
  }

  stats.data <- data.frame(Observed = init.statistics,
                           Simulated = colMeans(hsn.tot))
  rownames(stats.data) <- possible.stats
  cat("Statistics of observed network and networks simulated from final theta parameter estimates:\n")
  print(stats.data)

  statistic_test_p_values <- data.frame(statistic_test_p_values)
  rownames(statistic_test_p_values) <- possible.stats
  cat("\nt-test p values for statistics of observed network and networks simulated from final theta parameter estimates:\n \n")
  print(statistic_test_p_values)
  colnames(statistic_test_p_values) <- "p_values"

  #test to see if we have an acceptable fit
  acceptable_fit <- statistic_test_p_values[which(GERGM_Object@stats_to_use == 1),1]

  if(min(acceptable_fit) > acceptable_fit_p_value_threshold){
    GERGM_Object@acceptable_fit <- TRUE
    message("Parameter estimates simulate networks that are statistically indistinguishable from observed network. ")
  }else{
    GERGM_Object@acceptable_fit <- FALSE
    message("Parameter estimates simulate networks that are statistically distinguishable from observed network. Considder respecifying.")
  }

  # make GOF plot
  Gof_Plot(GERGM_Object)
  #4. Return GERGM object

  return(GERGM_Object)
}
