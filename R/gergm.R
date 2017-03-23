#' @title A Function to estimate a GERGM.
#' @description The main function provided by the package.
#'
#' @param formula A formula object that specifies the relationship between
#' statistics and the observed network. Currently, the user may specify a model
#' using any combination of the following statistics: `out2stars(alpha = 1)`,
#' `in2stars(alpha = 1)`, `ctriads(alpha = 1)`, `mutual(alpha = 1)`,
#' `ttriads(alpha = 1)`, `absdiff(covariate = "MyCov")`,
#'  `sender(covariate = "MyCov")`,
#' `reciever(covariate = "MyCov")`, `nodematch(covariate)`,
#' `nodemix(covariate, base = "MyBase")`, `netcov(network)` and
#' `edges(alpha = 1, method = c("regression","endogenous"))`. If the user
#' specifies `nodemix(covariate, base = NULL)`, then all levels of the
#' covariate will be matched on. Note that the `edges` term must be specified if
#' the user wishes to include an intercept
#' (strongly recommended). The user may select the "regression" method
#' (default) to include an intercept in the lambda transformation of the
#' network, or "endogenous" to include the intercept as in a traditional ERGM
#' model.  To use exponential down-weighting for any of the network level terms,
#' simply specify a value for alpha less than 1. The `(alpha = 1)` term may be
#' omitted from the structural terms if no exponential down weighting is
#' required. In this case, the terms may be provided as: `out2star`, `in2star`,
#' `ctriads`, `mutual`, `ttriads`. If the network is undirected the user may only
#' specify the following terms: `twostars(alpha = 1)`,  `ttriads(alpha = 1)`,
#' `absdiff(covariate = "MyCov")`,
#' `sender(covariate = "MyCov")`, `nodematch(covariate)`,
#' `nodemix(covariate, base = "MyBase")`, `netcov(network)` and
#' `edges(alpha = 1, method = c("regression","endogenous"))`. In some cases,
#' the user may only wish to calculate endogenous statistics for edges between
#' some subset of the nodes in the network. For each of the endogenous
#' statistics, the user may optionally specify a `covariate` and `base` field
#' such as in `in2stars(covariate = "Type", base = "C")`. This will add an
#' in-2star statistic for the subnetwork defined by actors who match each level
#' of the categorical variable "Type" (in this example), and exclude the
#' subnetwork for type "C", if the `base` argument is provided. If the `base`
#' argument is excluded, then terms will be added to the model for all levels of
#' the statistic. This can be a useful option if the user believes that a
#' network property varies with some property of nodes.
#' @param covariate_data A data frame containing node level covariates the user
#' wished to transform into sender or receiver effects. It must have row names
#' that match every entry in colnames(raw_network), should have descriptive
#' column names.  If left NULL, then no sender or receiver effects will be
#' added.
#' @param beta_correlation_model Defaults to FALSE. If TRUE, then the beta
#' correlation model is estimated. A correlation network must be provided, but
#' all covariates and undirected statistics may be supplied as normal.
#' @param distribution_estimator Provides an option to estimate the structure of
#' row-wise marginal and joint distribtuions using a uniform-dirichlet proposal
#' distribution. THIS FEATURE IS EXPERIMENTAL. Defaults to "none", in which case
#' a normal GERGM is estimated, but can be set to one of "rowwise-marginal" and
#' "joint" to propose either row-wsie marginal distribtuions or joint
#' distributions. If an option other than "none" is selected, the beta correlation
#' model will be turned off, estimation will automatically be set to Metropolis,
#' no covariate data will be allowed, and the network will be set to
#' directed. Furthermore, a "diagonal" statistic will be added to the model which
#' simply records the sum of the diagonal of the network. The "mutual" statistic
#' will also be adapted to include the diagonal elements. In the future, more
#' statistics which take account of the network diagonal will be included.
#' @param network_is_directed Logical specifying whether or not the observed
#' network is directed. Default is TRUE.
#' @param include_diagonal Logical indicating whether the diagonal should be
#' included in the estimation proceedure. If TRUE, then a "diagonal" statistic
#' is added to the model. Defaults to FALSE.
#' @param number_of_networks_to_simulate Number of simulations generated for
#' estimation via MCMC. Default is 500.
#' @param MCMC_burnin Number of samples from the MCMC simulation procedure that
#' will be discarded before drawing the samples used for estimation.
#' Default is 100.
#' @param thin The proportion of samples that are kept from each simulation. For
#' example, thin = 1/200 will keep every 200th network in the overall simulated
#' sample. Default is 1.
#' @param proposal_variance The variance specified for the Metropolis Hastings
#' simulation method. This parameter is inversely proportional to the average
#' acceptance rate of the M-H sampler and should be adjusted so that the average
#' acceptance rate is approximately 0.25. Default is 0.1.
#' @param target_accept_rate The target Metropolis Hastings acceptance rate.
#' Defaults to 0.25
#' @param seed Seed used for reproducibility. Default is 123.
#' @param hyperparameter_optimization Logical indicating whether automatic
#' hyperparameter optimization should be used. Defaults to FALSE. If TRUE, then
#' the algorithm will automatically seek to find an optimal burnin and number of
#' networks to simulate, and if using Metropolis Hasings, will attempt to select
#' a proposal variance that leads to a acceptance rate within +-0.05 of
#' target_accept_rate. Furthermore, if degeneracy is detected, the algorithm
#' will attempt to adress the issue automatically. WARNING: This feature is
#' experimental, and may greatly increase runtime. Please monitor console
#' output!
#' @param theta_grid_optimization_list Defaults to NULL. This highly
#' experimental feature may allow the user to address model degeneracy arising
#' from a suboptimal theta initialization. It performs a grid search around the
#' theta values calculated via MPLE to select a potentially improved
#' initialization. The runtime complexity of this feature grows exponentially in
#' the size of the grid and number of parameters -- use with great care. This
#' feature may only be used if hyperparameter_optimization = TRUE, and if a list
#' object of the following form is provided: list(grid_steps = 2,
#' step_size = 0.5, cores = 2, iteration_fraction = 0.5). grid_steps indicates
#' the number of steps out the grid search will perform, step_size indicates the
#' fraction of the MPLE theta estimate that each grid search step will change by,
#' cores indicates the number of cores to be used for parallel optimization, and
#' iteration_fraction indicates the fraction of the number of MCMC iterations
#' that will be used for each grid point (should be set less than 1 to speed up
#' optimization). In general grid_steps should be smaller the more structural
#' parameters the user wishes to specify. For example, with 5 structural
#' parameters (mutual, ttriads, etc.), grid_steps = 3 will result in a (2*3+1)^5
#' = 16807 parameter grid search. Again this feature is highly experimental and
#' should only be used as a last resort (after playing with exponential
#' down weighting and the MPLE_gain_factor).
#' @param weighted_MPLE Defaults to FALSE. Should be used whenever the user is
#' specifying statistics with alpha down weighting. Tends to provide better
#' initialization when downweight_statistics_together = FALSE.
#' @param parallel Logical indicating whether the weighted MPLE objective and any
#' other operations that can be easily parallelized should be calculated in
#' parallel. Defaults to FALSE. If TRUE, a significant speedup in computation
#' may be possible.
#' @param parallel_statistic_calculation Logical indicating whether network
#' statistics should be calculated in parallel. This will tend to be slower for
#' networks with les than ~30 nodes but may provide a substantial speedup for
#' larger networks.
#' @param cores Numeric value defaulting to 1. Can be set to any number up to the
#' number of threads/cores available on your machine. Will be used to speed up
#' computations if parallel = TRUE.
#' @param use_stochastic_MH A logical indicating whether a stochastic approximation
#' to the h statistics should be used under Metropolis Hastings in-between
#' thinned samples. This may dramatically speed up estimation. Defaults to FALSE.
#' HIGHLY EXPERIMENTAL!
#' @param stochastic_MH_proportion Percentage of dyads/triads to use for
#' approximation, defaults to 0.25.
#' @param slackr_integration_list An optional list object that contains
#' information necessary to provide updates about model fitting progress to a
#' Slack channel (https://slack.com/). This can be useful if models take a long
#' time to run, and you wish to receive updates on their progress (or if they
#' become degenerate). The list object must be of the following form:
#' list(model_name = "descriptive model name", channel = "#yourchannelname",
#'  incoming_webhook_url = "https://hooks.slack.com/services/XX/YY/ZZ"). You
#'  will need to set up incoming webhook integration for your slack channel and
#'  then paste in the URL you get from slack into the incoming_webhook_url field.
#'  If all goes well, and the computer you are running the GERGM estimation on
#'  has internet access, your slack channel will receive updates when you start
#'  estimation, after each lambda/theta parameter update, if the model becomes
#'  degenerate, and when it completes running.
#' @param convergence_tolerance Threshold designated for stopping criterion. If
#' the difference of parameter estimates from one iteration to the next all have
#' a p -value (under a paired t-test) greater than this value, the parameter
#' estimates are declared to have converged. Default is 0.5, which is quite
#' conservative.
#' @param MPLE_gain_factor Multiplicative constant between 0 and 1 that controls
#' how far away the initial theta estimates will be from the standard MPLEs via
#' a one step Fisher update. In the case of strongly dependent data, it is
#' suggested to use a value of 0.10. Default is 0.
#' @param acceptable_fit_p_value_threshold A p-value threshold for how closely
#' statistics of observed network conform to statistics of networks simulated
#' from GERGM parameterized by converged final parameter estimates. Default value
#' is 0.05.
#' @param normalization_type If only a raw_network is provided the function
#' will automatically check to determine if all edges fall in the [0,1] interval.
#' If edges are determined to fall outside of this interval, then a trasformation
#' onto the interval may be specified. If "division" is selected, then the data
#' will have a value added to them such that the minimum value is at least zero
#' (if necessary) and then all edge values will be divided by the maximum to
#' ensure that the maximum value is in [0,1]. If "log" is selected, then the data
#' will have a value added to them such that the minimum value is at least zero
#' (if necessary), then 1 will be added to all edge values before they are logged
#' and then divided by the largest value, again ensuring that the resulting
#' network is on [0,1]. Defaults to "log" and need not be set to NULL if
#' providing covariates as it will be ignored.
#' @param transformation_type Specifies how covariates are transformed onto the
#' raw network. When working with heavy tailed data that are not strictly
#' positive, select "Cauchy" to transform the data using a Cauchy distribution.
#' If data are strictly positive and heavy tailed (such as financial data) it is
#' suggested the user select "LogCauchy" to perform a Log-Cauchy transformation
#' of the data. For a tranformation of the data using a Gaussian distribution,
#' select "Gaussian" and for strictly positive raw networks, select "LogNormal".
#' The Default value is "Cauchy".
#' @param estimation_method Simulation method for MCMC estimation. Default is
#' "Metropolis", which allows for the most flexible model specifications, but
#' may also be set to "Gibbs", if the user wishes to use Gibbs sampling.
#' @param downweight_statistics_together Logical specifying whether or not the
#' weights should be applied inside or outside the sum. Default is TRUE and user
#' should not select FALSE under normal circumstances.
#' @param force_x_theta_updates Defaults to 1 where theta estimation is not
#' allowed to converge until thetas have updated for x iterations . Useful when
#' model is not degenerate but simulated statistics do not match observed network
#' well when algorithm stops after first y updates.
#' @param force_x_lambda_updates Defaults to 1 where lambda estimation is not
#' allowed to converge until lambdas have updated for x iterations . Useful when
#' model is not degenerate but simulated statistics do not match observed network
#' well when algorithm stops after first y updates.
#' @param output_directory The directory where you would like output generated
#' by the GERGM estimation procedure to be saved (if output_name is specified).
#' This includes, GOF, trace, and parameter estimate plots, as well as a summary
#' of the estimation procedure and an .Rdata file containing the GERGM object
#' returned by this function. May be left as NULL if the user would prefer all
#' plots be printed to the graphics device.
#' @param output_name The common name stem you would like to assign to all
#' objects output by the gergm() function. Default value of NULL will not save any
#' output directly to .pdf files, it will be printed to the console instead. Must
#' be a character string or NULL. For example, if "Test" is supplied as the
#' output_name, then 4 files will be output: "Test_GOF.pdf", "Test_Parameter_Estim
#' ates.pdf", "Test_GERGM_Object.Rdata", "Test_Estimation_Log.txt", and
#' "Test_Trace_Plot.pdf"
#' @param generate_plots Defaults to TRUE, if FALSE, then no diagnostic or
#' parameter plots are generated.
#' @param verbose Defaults to TRUE (providing lots of output while model is
#' running). Can be set to FALSE if the user wishes to see less output.
#' @param fine_grained_pv_optimization Logical indicating whether fine grained
#' proposal variance optimization should be used. This will often slow down
#' proposal variance optimization, but may provide better results. Highly
#' recommended if running a correlation model.
#' @param use_MPLE_only Logical specifying whether or not only the maximum pseudo
#' likelihood estimates should be obtained. In this case, no simulations will be
#' performed. Default is FALSE.
#' @param stop_for_degeneracy When TRUE, automatically stops estimation when
#' degeneracy is detected, even when hyperparameter_optimization is set to TRUE.
#' Defaults to FALSE.
#' @param maximum_number_of_lambda_updates Maximum number of iterations of outer
#' MCMC loop which alternately estimates transform parameters and ERGM
#' parameters. In the case that data_transformation = NULL, this argument is
#' ignored. Default is 10.
#' @param maximum_number_of_theta_updates Maximum number of iterations within the
#' MCMC inner loop which estimates the ERGM parameters. Default is 100.
#' @param estimate_model Logical indicating whether a model should be estimated.
#' Defaults to TRUE, but can be set to FALSE if the user simply wishes to return
#' a GERGM object containing the model specification. Useful for debugging.
#' @param ... Optional arguments, currently unsupported.
#' @return A gergm object containing parameter estimates.
#' @examples
#' \dontrun{
#' set.seed(12345)
#' net <- matrix(rnorm(100,0,20),10,10)
#' colnames(net) <- rownames(net) <- letters[1:10]
#' formula <- net ~  mutual + ttriads
#'
#' test <- gergm(formula,
#'               normalization_type = "division",
#'               network_is_directed = TRUE,
#'               use_MPLE_only = FALSE,
#'               estimation_method = "Metropolis",
#'               number_of_networks_to_simulate = 40000,
#'               thin = 1/10,
#'               proposal_variance = 0.5,
#'               downweight_statistics_together = TRUE,
#'               MCMC_burnin = 10000,
#'               seed = 456,
#'               convergence_tolerance = 0.01,
#'               MPLE_gain_factor = 0,
#'               force_x_theta_updates = 4)
#' }
#' @export
gergm <- function(formula,
                  covariate_data = NULL,
                  beta_correlation_model = FALSE,
                  distribution_estimator = c("none","rowwise-marginal","joint"),
                  network_is_directed = TRUE,
                  include_diagonal = FALSE,
                  number_of_networks_to_simulate = 500,
                  MCMC_burnin = 100,
                  thin = 1,
                  proposal_variance = 0.1,
                  target_accept_rate = 0.25,
                  seed = 123,
                  hyperparameter_optimization = FALSE,
                  theta_grid_optimization_list = NULL,
                  weighted_MPLE = FALSE,
                  parallel = FALSE,
                  parallel_statistic_calculation = FALSE,
                  cores = 1,
                  use_stochastic_MH = FALSE,
                  stochastic_MH_proportion = 0.25,
                  slackr_integration_list = NULL,
                  convergence_tolerance = 0.5,
                  MPLE_gain_factor = 0,
                  acceptable_fit_p_value_threshold = 0.05,
                  normalization_type = c("log","division"),
                  transformation_type = c("Cauchy","LogCauchy","Gaussian","LogNormal"),
                  estimation_method = c("Metropolis","Gibbs"),
                  downweight_statistics_together = TRUE,
                  force_x_theta_updates = 1,
                  force_x_lambda_updates = 1,
                  output_directory = NULL,
                  output_name = NULL,
                  generate_plots = TRUE,
                  verbose = TRUE,
                  fine_grained_pv_optimization = FALSE,
                  use_MPLE_only = c(FALSE, TRUE),
                  stop_for_degeneracy = FALSE,
                  maximum_number_of_lambda_updates = 10,
                  maximum_number_of_theta_updates = 10,
                  estimate_model = TRUE,
                  ...
                  ){

  # pass in experimental features through elipsis
  object <- as.list(substitute(list(...)))[-1L]

  # record start time for estimation
  start_time <- Sys.time()

  # hard coded possible stats
  possible_structural_terms <- c("out2stars",
                                 "in2stars",
                                 "ctriads",
                                 "mutual",
                                 "ttriads",
                                 "edges",
                                 "diagonal")
  possible_structural_terms_undirected <- c("twostars",
                                            "ttriads",
                                            "edges")
  possible_covariate_terms <- c("absdiff",
                                "nodecov",
                                "nodematch",
                                "sender",
                                "receiver",
                                "intercept",
                                "nodemix")
  possible_network_terms <- "netcov"
  possible_transformations <- c("cauchy",
                                "logcauchy",
                                "gaussian",
                                "lognormal")


  # set logical values for whether we are using MPLE only, whether the network
  # is directed, and which estimation method we are using as well as the
  # transformation type
  use_MPLE_only <- use_MPLE_only[1] #default is FALSE
  network_is_directed <- network_is_directed[1] #default is TRUE
  estimation_method <- estimation_method[1] #default is MH
  transformation_type <- transformation_type[1] #default is "Cauchy"
  transformation_type <- tolower(transformation_type)
  normalization_type <- normalization_type[1]
  distribution_estimator <- distribution_estimator[1]
  using_distribution_estimator <- FALSE

  # deal with the case where we are using a distribution estimator
  if (distribution_estimator %in%  c("none","rowwise-marginal","joint")) {
    # if we are actually using the distribution estimator
    if (distribution_estimator != "none") {
      # perform checks and set variables so that they are at their correct values.
      cat("Making sure all options are set properly for use with the",
          "distribution estimator... \n")
      use_MPLE_only <- FALSE
      estimation_method <- "Metropolis"
      covariate_data <- NULL
      network_is_directed <- TRUE
      normalization_type <- "division"
      beta_correlation_model <- FALSE
      using_distribution_estimator <- TRUE
      include_diagonal <- TRUE
    }
  } else {
    stop("distribution_estimator must be one of 'none','rowwise-marginal', or 'joint'")
  }

  # add a diagonal term if requested.
  if (include_diagonal) {
    if (estimation_method == "Gibbs") {
      cat("Gibbs sampling is currently not supported when include_diagonal == TRUE, changing to 'Metropolis'")
      estimation_method <- "Metropolis"
    }
    formula <- add_diagonal_term(formula)
  }

  # set the number of threads to use with parallel
  if (parallel) {
    RcppParallel::setThreadOptions(numThreads = cores)
  }

  # if we are using a correlation network, then the network must be undirected.
  if (beta_correlation_model) {
    if (network_is_directed) {
      cat("Setting network_is_directed to FALSE for correlation network...\n")
    }
    network_is_directed <- FALSE
  }

  if (network_is_directed) {
    possible_structural_term_indices <- 1:7
  } else {
    possible_structural_term_indices <- c(2,5,6,7)
  }

  # check terms for undirected network
  if (!network_is_directed) {
    formula <- parse_undirected_structural_terms(
      formula,
      possible_structural_terms,
      possible_structural_terms_undirected)
  }

  # check to see if we are using a regression intercept term, and if we are,
  # then add the "intercept" field to the formula.
  check_for_edges <- Parse_Formula_Object(formula,
                               possible_structural_terms,
                               possible_covariate_terms,
                               possible_network_terms,
                               using_distribution_estimator,
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

  if (is.null(output_directory) & !is.null(output_name)) {
    stop("You have specified an output file name but no output directory. Please
         specify both or neither.")
  }

  if (length(which(possible_transformations %in% transformation_type  == T)) != 1) {
    stop("You have specified a transformation that is not recognized. Please
         specify one of: Cauchy, LogCauchy, Gaussian, or LogNormal")
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
     is_correlation_network = FALSE,
     is_directed = network_is_directed,
     beta_correlation_model = beta_correlation_model,
     include_diagonal = include_diagonal)

  data_transformation <- NULL
  if (!is.null(Transformed_Data$transformed_covariates)) {
    data_transformation <- Transformed_Data$transformed_covariates
  }
  gpar.names <- c(Transformed_Data$gpar.names, "dispersion")

  #1. Create GERGM object from network

  GERGM_Object <- Create_GERGM_Object_From_Formula(
     formula,
     theta.coef = NULL,
     possible_structural_terms,
     possible_covariate_terms,
     possible_network_terms,
     using_distribution_estimator,
     raw_network = Transformed_Data$network,
     together = 1,
     transform.data = data_transformation,
     lambda.coef = NULL,
     transformation_type = transformation_type,
     is_correlation_network = FALSE,
     is_directed = network_is_directed,
     beta_correlation_model = beta_correlation_model,
     covariate_data = covariate_data,
     possible_structural_terms_undirected = possible_structural_terms_undirected)

  GERGM_Object@theta_estimation_converged <- FALSE
  GERGM_Object@lambda_estimation_converged <- FALSE
  GERGM_Object@observed_network  <- GERGM_Object@network
  GERGM_Object@observed_bounded_network <- GERGM_Object@bounded.network
  GERGM_Object@simulation_only <- FALSE
  GERGM_Object@transformation_type <- transformation_type
  GERGM_Object@downweight_statistics_together <- downweight_statistics_together
  GERGM_Object@directed_network <- network_is_directed
  # only add in list if not NULL
  GERGM_Object@using_grid_optimization <- FALSE
  if (class(theta_grid_optimization_list) == "list") {
    GERGM_Object@using_grid_optimization <- TRUE
    GERGM_Object@theta_grid_optimization_list <- theta_grid_optimization_list
  }

  if (!is.null(data_transformation)) {
    GERGM_Object@data_transformation <- data_transformation
  }

  if (is.null(output_name)) {
    GERGM_Object@print_output <- FALSE
  }else{
    GERGM_Object@print_output <- TRUE
  }

  # if we are using a correlation network then set field to TRUE.
  GERGM_Object@is_correlation_network <- FALSE # deprecated
  GERGM_Object@beta_correlation_model <- beta_correlation_model
  GERGM_Object@distribution_estimator <- distribution_estimator
  GERGM_Object@include_diagonal <- include_diagonal

  # record the various optimizations we are using so that they can be used in
  # the main algorithm
  GERGM_Object@weighted_MPLE <- weighted_MPLE
  GERGM_Object@fine_grained_pv_optimization <- fine_grained_pv_optimization
  GERGM_Object@parallel <- parallel
  GERGM_Object@parallel_statistic_calculation <- parallel_statistic_calculation
  GERGM_Object@cores <- cores
  GERGM_Object@use_stochastic_MH <- use_stochastic_MH
  GERGM_Object@stochastic_MH_proportion <- stochastic_MH_proportion
  GERGM_Object@possible_endogenous_statistic_indices <- possible_structural_term_indices

  # set adaptive metropolis parameters
  GERGM_Object@hyperparameter_optimization <- hyperparameter_optimization
  GERGM_Object@target_accept_rate <- target_accept_rate
  GERGM_Object@proposal_variance <- proposal_variance
  GERGM_Object@estimation_method <- estimation_method
  GERGM_Object@number_of_simulations <- number_of_networks_to_simulate
  GERGM_Object@thin <- thin
  GERGM_Object@burnin <- MCMC_burnin
  GERGM_Object@MPLE_gain_factor <- MPLE_gain_factor
  GERGM_Object@start_time <- toString(start_time)

  GERGM_Object@using_slackr_integration <- FALSE
  if (!is.null(slackr_integration_list)) {
    if (length(slackr_integration_list) != 3) {
      stop("slackr_integration_list must contain three fields, see documentation...")
    }
    GERGM_Object@slackr_integration_list <- slackr_integration_list
    GERGM_Object@using_slackr_integration <- TRUE

    start_message <- paste('Starting GERGM Estimation at',toString(start_time))
    specification <- paste('Specification:',toString(formula))
    specification <- gsub("\"", "", specification, fixed=TRUE)

    slackr::slackr_bot(
      start_message,
      specification,
      channel = GERGM_Object@slackr_integration_list$channel,
      username = GERGM_Object@slackr_integration_list$model_name,
      incoming_webhook_url = GERGM_Object@slackr_integration_list$incoming_webhook_url)
  }

  # prepare auxiliary data
  GERGM_Object@statistic_auxiliary_data <- prepare_statistic_auxiliary_data(
    GERGM_Object)

  #GERGM_Object@statistic_auxiliary_data$triples
  #GERGM_Object@statistic_auxiliary_data$pairs
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

  #2. Estimate GERGM
  if (estimate_model) {
    GERGM_Object <- Estimate_GERGM(
      formula,
      MPLE.only = use_MPLE_only,
      max.num.iterations = maximum_number_of_lambda_updates,
      mc.num.iterations = maximum_number_of_theta_updates,
      seed = seed,
      tolerance = convergence_tolerance,
      possible.stats = possible_structural_terms,
      GERGM_Object = GERGM_Object,
      force_x_theta_updates = force_x_theta_updates,
      transformation_type =  transformation_type,
      verbose = verbose,
      force_x_lambda_updates = force_x_lambda_updates,
      stop_for_degeneracy = stop_for_degeneracy)

    #3. Perform degeneracy diagnostics and create GOF plots
    if (!GERGM_Object@theta_estimation_converged) {
      warning("Estimation procedure did not detect convergence in Theta estimates. Estimation halted when maximum number of updates was reached. Be careful to assure good model fit or select a more relaxed convergence criterion.")
      GERGM_Object <- store_console_output(GERGM_Object,"Estimation procedure did not detect convergence in Theta estimates. Estimation halted when maximum number of updates was reached. Be careful to assure good model fit or select a more relaxed convergence criterion.")
    }
    if (!GERGM_Object@lambda_estimation_converged) {
      warning("Estimation procedure did not detect convergence in Lambda estimates. Estimation halted when maximum number of updates was reached. Be careful to assure good model fit or select a more relaxed convergence criterion.")
      GERGM_Object <- store_console_output(GERGM_Object,"Estimation procedure did not detect convergence in Lambda estimates. Estimation halted when maximum number of updates was reached. Be careful to assure good model fit or select a more relaxed convergence criterion.")
    }

    #now simulate from last update of theta parameters
    GERGM_Object <- Simulate_GERGM(
      GERGM_Object,
      seed1 = seed,
      possible.stats = possible_structural_terms,
      parallel = GERGM_Object@parallel_statistic_calculation)

    colnames(GERGM_Object@lambda.coef) = gpar.names

    # change back column names if we are dealing with an undirected network
    if (!network_is_directed) {
      change <- which(colnames(GERGM_Object@theta.coef) == "in2stars")
      if (length(change) > 0) {
        colnames(GERGM_Object@theta.coef)[change] <- "twostars"
      }
    }

    init.statistics <- NULL
    if (GERGM_Object@is_correlation_network) {
      init.statistics <- calculate_h_statistics(
        GERGM_Object,
        GERGM_Object@statistic_auxiliary_data,
        all_weights_are_one = FALSE,
        calculate_all_statistics = TRUE,
        use_constrained_network = FALSE)
    }else{
      init.statistics <- calculate_h_statistics(
        GERGM_Object,
        GERGM_Object@statistic_auxiliary_data,
        all_weights_are_one = FALSE,
        calculate_all_statistics = TRUE,
        use_constrained_network = TRUE)
    }
    # fix issue with the wrong stats being saved
    GERGM_Object@stats[2,] <- init.statistics
    hsn.tot <- GERGM_Object@MCMC_output$Statistics

    # save these statistics so we can make GOF plots in the future, otherwise
    # they would be the transformed statistics which would produce poor GOF plots.
    GERGM_Object@simulated_statistics_for_GOF <- hsn.tot
    #thin statsitics
    hsn.tot <- Thin_Statistic_Samples(hsn.tot)

    #calculate t.test p-values for calculating the difference in the means of
    # the newly simulated data with the original network
    statistic_test_p_values <- rep(NA,ncol(hsn.tot))
    for (i in 1:ncol(hsn.tot)) {
      if (length(unique(hsn.tot[,i])) > 1) {
        statistic_test_p_values[i] <- round(t.test(hsn.tot[, i],
                                                   mu = init.statistics[i])$p.value,3)
      }
    }


    stats.data <- data.frame(Observed = init.statistics,
                             Simulated = colMeans(hsn.tot))
    rownames(stats.data) <- GERGM_Object@full_theta_names
    cat("Statistics of observed network and networks simulated from final theta parameter estimates:\n")
    GERGM_Object <- store_console_output(GERGM_Object,"Statistics of observed network and networks simulated from final theta parameter estimates:\n")

    GERGM_Object <- store_console_output(GERGM_Object, toString(stats.data))

    statistic_test_p_values <- data.frame(p_values = statistic_test_p_values)
    rownames(statistic_test_p_values) <- GERGM_Object@full_theta_names
    cat("\nt-test p-values for statistics of observed network and networks simulated from final theta parameter estimates:\n \n")
    GERGM_Object <- store_console_output(GERGM_Object,"\nt-test p-values for statistics of observed network and networks simulated from final theta parameter estimates:\n \n")
    print(statistic_test_p_values)
    GERGM_Object <- store_console_output(GERGM_Object, toString(statistic_test_p_values))


    colnames(statistic_test_p_values) <- "p_values"
    GERGM_Object@observed_simulated_t_test <- statistic_test_p_values

    #test to see if we have an acceptable fit
    check_stats <- GERGM_Object@statistic_auxiliary_data$specified_statistic_indexes_in_full_statistics
    acceptable_fit <- statistic_test_p_values[check_stats, 1]

    if (min(acceptable_fit) > acceptable_fit_p_value_threshold) {
      GERGM_Object@acceptable_fit <- TRUE
      message("Parameter estimates simulate networks that are statistically indistinguishable from observed network on the statistics specified by the user. ")
      GERGM_Object <- store_console_output(GERGM_Object,"Parameter estimates simulate networks that are statistically indistinguishable from observed network on the statistics specified by the user. ")
    }else{
      GERGM_Object@acceptable_fit <- FALSE
      message("Parameter estimates simulate networks that are statistically distinguishable from observed network. Check GOF plots to determine if the model provides a reasonable fit . This is a very stringent test for goodness of fit, so results may still be acceptable even if this criterion is not met.")
      GERGM_Object <- store_console_output(GERGM_Object, "Parameter estimates simulate networks that are statistically distinguishable from observed network. Check GOF plots to determine if the model provides a reasonable fit . This is a very stringent test for goodness of fit, so results may still be acceptable even if this criterion is not met.")
    }

    GERGM_Object@simulated_bounded_networks_for_GOF <- GERGM_Object@MCMC_output$Networks

    # precalculate intensities and degree distributions so we can return them.
    GERGM_Object <- calculate_additional_GOF_statistics(GERGM_Object)

    #4. output everything to the appropriate files and return GERGM object.
    if (generate_plots) {
      # only generate output if output_name is not NULL
      if (!is.null(output_name)) {
        if (is.null(output_directory)) {
          output_directory <- getwd()
        }
        current_directory <- getwd()
        setwd(output_directory)

        try({
          pdf(file = paste(output_name,"_GOF.pdf",sep = ""),
              height = 4,
              width = 10)
          GOF(GERGM_Object)
          dev.off()

          pdf(file = paste(output_name,"_Parameter_Estimates.pdf",sep = ""),
              height = 4,
              width = 5)
          Estimate_Plot(GERGM_Object)
          dev.off()

          pdf(file = paste(output_name,"_Trace_Plot.pdf",sep = ""),
              height = 4,
              width = 6)
          Trace_Plot(GERGM_Object)
          dev.off()

          save(GERGM_Object,
               file = paste(output_name,"_GERGM_Object.RData",sep = ""))

          write.table(GERGM_Object@console_output,
                      file = paste(output_name,"_Estimation_Log.txt",sep = ""),
                      row.names = F,
                      col.names = F,
                      fileEncoding = "utf8",
                      quote = F)
        })
        setwd(current_directory)
      } else{
        # if we are not saving everything to a directory then just print stuff to
        # the graphics device
        try({
          GOF(GERGM_Object)
          Sys.sleep(2)
          Estimate_Plot(GERGM_Object)
          Sys.sleep(2)
          Trace_Plot(GERGM_Object)
        })
      }
    }

    # transform networks back to observed scale
    cat("Transforming networks simulated via MCMC as part of the fit diagnostics back on to the scale of observed network. You can access these networks through the '@MCMC_output$Networks' field returned by this function...\n")
    GERGM_Object <- Convert_Simulated_Networks_To_Observed_Scale(GERGM_Object)
  } else {
    cat("Returning GERGM object without estimating model...\n")
  }

  # now report the total elapsed time and end time to the user
  end_time <- Sys.time()
  GERGM_Object@end_time <- toString(end_time)

  elapsed_time <- end_time - start_time
  GERGM_Object@elapsed_time <- toString(elapsed_time)

  cat("Estimation Complete at:",toString(end_time),
      "\nElapsed time (seconds):",elapsed_time,"\n")

  # push to slack if desired
  if (GERGM_Object@using_slackr_integration) {
    completion_time <- paste("Estimation Complete at:",toString(end_time))
    elapsed_time <- paste("Elapsed time (seconds):",elapsed_time)
    slackr::slackr_bot(
      completion_time,
      elapsed_time,
      channel = GERGM_Object@slackr_integration_list$channel,
      username = GERGM_Object@slackr_integration_list$model_name,
      incoming_webhook_url = GERGM_Object@slackr_integration_list$incoming_webhook_url)
  }

  return(GERGM_Object)
}
