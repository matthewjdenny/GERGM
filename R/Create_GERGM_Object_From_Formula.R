# create a gergm from a formula object (getgergm)
Create_GERGM_Object_From_Formula <- function(object,
                                             theta.coef,
                                             possible_structural_terms,
                                             possible_covariate_terms,
                                             possible_network_terms,
                                             using_distribution_estimator,
                                             raw_network,
                                             together = 1,
                                             transform.data = NULL,
                                             lambda.coef = NULL,
                                             transformation_type,
                                             is_correlation_network = FALSE,
                                             is_directed = TRUE,
                                             beta_correlation_model = FALSE,
                                             covariate_data = NULL,
                                             possible_structural_terms_undirected = NULL
                                             ){

  res1 <- Parse_Formula_Object(object,
                               possible_structural_terms,
                               possible_covariate_terms,
                               possible_network_terms,
                               using_distribution_estimator = using_distribution_estimator,
                               raw_network = raw_network,
                               theta = theta.coef,
                               terms_to_parse = "structural",
                               covariate_data = covariate_data)
  thetas <- res1$thetas
  network <- res1$net
  alphas <- res1$alphas
  statistics <- res1$statistics
  thresholds <- res1$thresholds
  nodes <- res1$nodes
  non_base_statistic <- res1$non_base_statistic
  legacy_statistics <- res1$legacy_statistics
  legacy_alphas <- res1$legacy_alphas
  legacy_thetas <- res1$legacy_thetas
  legacy_thresholds <- res1$legacy_thresholds
  theta_names <- res1$theta_names
  covariate_terms_only <- res1$covariate_terms_only

  # get the full names so we can name the h.statistics correctly after simulation
  full_theta_names <- possible_structural_terms
  if (!is_directed) {
    full_theta_names <- possible_structural_terms_undirected
  }
  nbinds <- which(non_base_statistic == 1)
  if (length(nbinds) > 0) {
    full_theta_names <- c(full_theta_names,theta_names[nbinds])
  }

  # for now we are not going to allow any covariates
  if (is_correlation_network) {
    if (!is.null(lambda.coef)) {
      stop("Covariate effects are currently not supported for correlation networks. Please respecify without covariates.")
    }
  } else if (beta_correlation_model) {
    cat("Using Beta model for correlation network data...\n")
    # if we are using the beta correlation model
    diag(network) <- 1
    bounded.network <- correlations.to.partials(network)
  } else if (!is.null(lambda.coef)) {
    cat("Covariates Provided...\n")
    # create the network based on the transform family
    # if there are no lambda.coefficients, we assume there is no transformation
    # if there is a transformation specified, transform the observed network
    if (transformation_type == "logcauchy" | transformation_type == "lognormal") {
      if (min(network) <= 0) {
        stop(paste("You have selected either a log-Cauchy or log-normal transformation but you have provided a network with values that are less than or equal to zero. Please ensure that the minimum value of the network you provide is greater than zero, or select a cauchy or normal transformation. The minimum value of the network provided is:",min(network)))
      }
      network <- log(network)
    }
    beta <- lambda.coef[1:(length(lambda.coef) - 1)]
    sig <- 0.01 + exp(lambda.coef[length(lambda.coef)])
    BZ <- 0
    if (is.na(dim(transform.data)[3])) {
      BZ = BZ + beta * transform.data
    }
    if (!is.na(dim(transform.data)[3])) {
      for (j in 1:(dim(transform.data)[3])) {
        BZ <- BZ + beta[j] * transform.data[, , j]
      }
    }
    if (transformation_type == "logcauchy" | transformation_type == "cauchy") {
      bounded.network <- pst(network, BZ, sig, 1)
    }
    if (transformation_type == "lognormal" | transformation_type == "gaussian") {
      bounded.network <- pst(network, BZ, sig, Inf)
    }

  } # end of lambda conditional

  if (is.null(lambda.coef)) {
    bounded.network <- network
    lambda.coef <- as.data.frame(0)
  }

  if (!is.null(lambda.coef)) {
    lambda.coef <- as.data.frame(rbind(lambda.coef,NA))
    rownames(lambda.coef) <- c("est", "se")
  }

  # if we are providing a correlation network, transform it
  if (is_correlation_network) {
    diag(network) <- 1
    print(round(network,2))
    bounded.network <- transform.correlations(network)
  }

  theta.par <- thetas
  thetas <- t(as.matrix(thetas))
  thetas <- rbind(thetas, NA)
  colnames(thetas) <- theta_names
  rownames(thetas) <- c("est", "se")
  thetas <- as.data.frame(thetas)

  legacy_thetas <- t(as.matrix(legacy_thetas))
  legacy_thetas <- rbind(legacy_thetas, NA)
  colnames(legacy_thetas) <- possible_structural_terms
  rownames(legacy_thetas) <- c("est", "se")
  legacy_thetas <- as.data.frame(legacy_thetas)

  object <- Create_GERGM_Object(network = network,
                                bounded.network = bounded.network,
                                formula = object,
                                thetas = thetas,
                                lambda = lambda.coef,
                                alpha = alphas,
                                together = together,
                                possible.stats = possible_structural_terms,
                                thresholds = thresholds)

  object@stats_to_use <- statistics
  object@endogenous_statistic_node_sets <- nodes
  object@non_base_statistic_indicator  <- non_base_statistic
  object@legacy_statistics <- legacy_statistics
  object@legacy_alphas <-  legacy_alphas
  object@legacy_thetas <- legacy_thetas
  object@legacy_thresholds <- legacy_thresholds
  object@theta_names <- theta_names
  object@theta.par <- theta.par
  object@full_theta_names <- full_theta_names
  object@covariate_terms_only <- covariate_terms_only
  return(object)
}
