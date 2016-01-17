#' @title A Function to estimate a number of GERGMs in parallel, each with its own equation.
#' @description Allows the user to run multiple specifications at once in
#' parallel. All varaibles (excluding formula_list, observed_network_list,
#' covariate_data_list, network_data_list and cores) be be either specified as a
#' single value or as a vector of values equal to the length of formula_list, if
#' the user wishes to use different values for each specification. Currently
#' only functional for Unix based operating systems (OSX, Linux, CentOS, etc.)
#'
#' @param formula_list A list of formula objects that specifies the relationship
#' between statistics and the observed network for each gergm. See the gergm()
#' documentation for more details.
#' @param observed_network_list A list of observed networks (as numeric matrices
#' to be used with each specification).
#' @param covariate_data_list An optional list of covariate data frames (may
#' include NULL entries if no covariates are needed in some specifications)
#' @param network_data_list An optional list of of lists of network covariates
#' to be included in each specification (one list per specification -- may also be left NULL for some specifications). THe list object corresponding to each
#' specification must have entries for network covariates named as they appear
#' in the corresponding equation. For example if the user specified a
#' 'netcov(distance)' term, the corresponding list object for that specification
#' would need a $distance entry containing the corresponding matrix object.
#' @param cores The number of cores to be used for parallelization.
#' @param normalization_type If only a raw_network is provided and
#' omit_intercept_term = TRUE then, the function
#' will automatically check to determine if all edges fall in the [0,1] interval.
#' If edges are determined to fall outside of this interval, then a trasformation
#' onto the interval may be specified. If "division" is selected, then the data
#' will have a value added to them such that the minimum value is atleast zero
#' (if necessary) and then all edge values will be divided by the maximum to
#' ensure that the maximum value is in [0,1]. If "log" is selected, then the data
#' will have a value added to them such that the minimum value is atleast zero
#' (if necessary), then 1 will be added to all edge values before they are logged
#' and then divided by the largest value, again ensuring that the resulting
#' network is on [0,1]. Defaults to "log" and need not be set to NULL if
#' providing covariates as it will be ignored.
#' @param network_is_directed Logical specifying whether or not the observed
#' network is directed. Default is TRUE.
#' @param use_MPLE_only Logical specifying whether or not only the maximum pseudo
#' likelihood estimates should be obtained. In this case, no simulations will be
#' performed. Default is FALSE.
#' @param transformation_type Specifies how covariates are transformed onto the
#' raw network. When working with heavly tailed data that are not strictly
#' positive, select "Cauchy" to transform the data using a Cauchy distribution.
#' If data are strictly positive and heavy tailed (such as financial data) it is
#' suggested the user select "LogCauchy" to perform a Log-Cauchy transformation
#' of the data. For a tranformation of the data using a Gaussian distribution,
#' select "Gaussian" and for strictly positive raw networks, select "LogNormal".
#' The Default value is "Cauchy".
#' @param estimation_method Simulation method for MCMC estimation. Default is
#' "Gibbs" which will generally be faster with well behaved networks but will not
#' allow for exponential downweighting.
#' @param maximum_number_of_lambda_updates Maximum number of iterations of outer
#' MCMC loop which alternately estimates transform parameters and ERGM
#' parameters. In the case that data_transformation = NULL, this argument is
#' ignored. Default is 10.
#' @param maximum_number_of_theta_updates Maximum number of iterations within the
#' MCMC inner loop which estimates the ERGM parameters. Default is 100.
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
#' will be discarded before drawing the samples used for estimation.
#' Default is 100.
#' @param seed Seed used for reproducibility. Default is 123.
#' @param convergence_tolerance Threshold designated for stopping criterion. If
#' the difference of parameter estimates from one iteration to the next all have
#' a p -value (under a paired t-test) greater than this value, the parameter
#' estimates are declared to have converged. Default is 0.01.
#' @param MPLE_gain_factor Multiplicative constant between 0 and 1 that controls
#' how far away the initial theta estimates will be from the standard MPLEs via
#' a one step Fisher update. In the case of strongly dependent data, it is
#' suggested to use a value of 0.10. Default is 0.
#' @param acceptable_fit_p_value_threshold A p-value threshold for how closely
#' statistics of observed network conform to statistics of networks simulated
#' from GERGM parameterized by converged final parameter estimates. Default value
#' is 0.05.
#' @param force_x_theta_updates Defaults to 1 where theta estimation is not
#' allowed to converge until thetas have updated for x iterations . Useful when
#' model is not degenerate but simulated statistics do not match observed network
#' well when algorithm stops after first y updates.
#' @param force_x_lambda_updates Defaults to 1 where lambda estimation is not
#' allowed to converge until lambdas have updated for x iterations . Useful when
#' model is not degenerate but simulated statistics do not match observed network
#' well when algorithm stops after first y updates.
#' @param output_directory The directory where you would like output generated
#' by the GERGM estimation proceedure to be saved (if output_name is specified).
#' This includes, GOF, trace, and parameter estimate plots, as well as a summary
#' of the estimation proceedure and an .Rdata file containing the GERGM object
#' returned by this function. May be left as NULL if the user would prefer all
#' plots be printed to the graphics device.
#' @param output_name The common name stem you would like to assign to all
#' objects output by the gergm function. Default value of NULL will not save any
#' output directly to .pdf files, it will be printed to the console instead. Must
#' be a character string or NULL. For example, if "Test" is supplied as the
#' output_name, then 4 files will be output: "Test_GOF.pdf", "Test_Parameter_Estim
#' ates.pdf", "Test_GERGM_Object.Rdata", "Test_Estimation_Log.txt", and
#' "Test_Trace_Plot.pdf"
#' @param generate_plots Defaults to TRUE, if FALSE, then no diagnostic or
#' parameter plots are generated.
#' @param verbose Defaults to TRUE (providing lots of output while model is
#' running). Can be set to FALSE if the user wishes to see less output.
#' @param omit_intercept_term Defualts to FALSE, can be set to TRUE if the
#' user wishes to omit the model intercept term.
#' @param hyperparameter_optimization Logical indicating whether automatic
#' hyperparameter optimization should be used. Defaults to FALSE. If TRUE, then
#' the algorithm will automatically seek to find an optimal burnin and number of
#' networks to simulate, and if using Metropolis Hasings, will attempt to select
#' a proposal variance that leads to a acceptance rate within +-0.05 of
#' target_accept_rate. Furthermore, if degeneracy is detected, the algorithm
#' will attempt to adress the issue automatically. WARNING: This feature is
#' experimental, and may greatly increase runtime. Please monitor console
#' output!
#' @param target_accept_rate The target Metropolis Hastings acceptance rate.
#' Defaults to 0.25
#' @param ... Optional arguments, currently unsupported.
#' @return A list of gergm objects for each model specified.
#' @export
parallel_gergm <- function(
  formula_list,
  observed_network_list,
  covariate_data_list = NULL,
  network_data_list = NULL,
  cores = 1,
  normalization_type = c("log","division"),
  network_is_directed = c(TRUE, FALSE),
  use_MPLE_only = c(FALSE, TRUE),
  transformation_type = c("Cauchy","LogCauchy","Gaussian","LogNormal"),
  estimation_method = c("Gibbs", "Metropolis"),
  maximum_number_of_lambda_updates = 10,
  maximum_number_of_theta_updates = 10,
  number_of_networks_to_simulate = 500,
  thin = 1,
  proposal_variance = 0.1,
  downweight_statistics_together = TRUE,
  MCMC_burnin = 100,
  seed = 123,
  convergence_tolerance = 0.01,
  MPLE_gain_factor = 0,
  acceptable_fit_p_value_threshold = 0.05,
  force_x_theta_updates = 1,
  force_x_lambda_updates = 1,
  output_directory = NULL,
  output_name = NULL,
  generate_plots = TRUE,
  verbose = TRUE,
  omit_intercept_term = FALSE,
  hyperparameter_optimization = FALSE,
  target_accept_rate = 0.25,
  ...
){

  # get the number of specifications
  num_specifications <- length(formula_list)

  if(length(formula_list) == 1){
    stop("You have only included one specification, use the gergm() function.")
  }

  vec <- 1:num_specifications
  cat("Running",num_specifications,"GERGM specifications on",cores,
      "cores. This may take a while...\n")
  GERGM_Results_List <- parallel::mclapply(X = vec,
    FUN = hysteresis_parallel,
    mc.cores = cores,
    num_specifications = num_specifications,
    formula_list = formula_list,
    observed_network_list = observed_network_list,
    covariate_data_list = covariate_data_list,
    network_data_list = network_data_list,
    normalization_type = normalization_type,
    network_is_directed = network_is_directed,
    use_MPLE_only = use_MPLE_only,
    transformation_type = transformation_type,
    estimation_method = estimation_method,
    maximum_number_of_lambda_updates = maximum_number_of_lambda_updates,
    maximum_number_of_theta_updates = maximum_number_of_theta_updates,
    number_of_networks_to_simulate = number_of_networks_to_simulate,
    thin = thin,
    proposal_variance = proposal_variance,
    downweight_statistics_together = downweight_statistics_together,
    MCMC_burnin = MCMC_burnin,
    seed = seed,
    convergence_tolerance = convergence_tolerance,
    MPLE_gain_factor = MPLE_gain_factor,
    acceptable_fit_p_value_threshold = acceptable_fit_p_value_threshold,
    force_x_theta_updates = force_x_theta_updates,
    force_x_lambda_updates = force_x_lambda_updates,
    output_directory = output_directory,
    output_name = output_name,
    generate_plots = generate_plots,
    verbose = verbose,
    omit_intercept_term = omit_intercept_term,
    hyperparameter_optimization = hyperparameter_optimization,
    target_accept_rate = target_accept_rate,
    ... = ...)


  return(GERGM_Results_List)
  }
