Parse_Formula_Object <- function(formula,
                                 possible_structural_terms,
                                 possible_covariate_terms,
                                 possible_network_terms,
                                 using_distribution_estimator = FALSE,
                                 raw_network = NULL,
                                 theta = NULL,
                                 terms_to_parse = "structural",
                                 covariate_data = NULL) {

  num_stats <- length(possible_structural_terms)
  # parse the formula
  if (class(formula) != "formula") {
    stop("'formula' must be a formula object.")
  }
  lhs <- deparse(formula[[2]])  # name of the response variable

  # get the actual response data
  if(is.null(raw_network)){
#     temp_net <- mget(lhs, envir = environment(),ifnotfound = list(not_found = NA))
#     if(class(temp_net) == "list"){
#       temp_net <- mget(lhs, envir = globalenv(),ifnotfound = list(not_found = NA))
#     }
    net <- dynGet(as.character(lhs),
                  ifnotfound = get(as.character(lhs)))
  }else{
    net <- raw_network
  }

  # check that the response is a matrix
  if (class(net)[1] != "matrix") {
    stop("Response must be a matrix object.")
  }

  rhs <- paste0(deparse(formula[[3]]), collapse = "")  # rhs of formula
  rhs <- gsub("\\s+", "", rhs)  # get rid of redundant spaces
  rhs <- strsplit(rhs, "\\+")[[1]]  # parse separate formula elements
  parsed_rhs <- vector (length = length(rhs), mode = "list")
  rhs_term_names <- rep("", length(rhs))
  alpha <- rep(1, length(rhs))
  threshold <- rep(0, length(rhs))
  for (i in 1:length(rhs)){
    parsed_rhs[[i]] <- parse_formula_term(rhs[i],
                                         possible_structural_terms,
                                         possible_covariate_terms,
                                         possible_network_terms)
    rhs_term_names[i] <- parsed_rhs[[i]]$term
    alpha[i] <- as.numeric(parsed_rhs[[i]]$weight)
    threshold[i] <- as.numeric(parsed_rhs[[i]]$threshold)
  }

  # check to see if the edges method is equal to "endogenous", and if not,
  # remove it.
  remove_ind <- NULL
  regression_intercept <- FALSE
  include_intercept <- FALSE
  for (i in 1:length(rhs)){
    if (rhs_term_names[i] == "edges") {
      include_intercept <- TRUE
      # if we are using the distribtuion estimator, then make sure that we set
      # the edges term to endogenous
      if (using_distribution_estimator) {
        if (parsed_rhs[[i]]$method != "endogenous") {
          parsed_rhs[[i]]$method <- "endogenous"
        }
      } else {
        if (parsed_rhs[[i]]$method != "endogenous") {
          remove_ind <- i
          regression_intercept <- TRUE
        }
      }
    }
  }

  if(!is.null(remove_ind)) {
    parsed_rhs <- parsed_rhs[-remove_ind]
    rhs_term_names <- rhs_term_names[-remove_ind]
    alpha <- alpha[-remove_ind]
    threshold <- threshold[-remove_ind]
  }


  # if we are parsing the structural terms out of the formula
  if(terms_to_parse == "structural"){
    # remove all node level covariate terms

    # deal with the case where we have zero structural covariates
    if (length(parsed_rhs) > 0) {
      remove <- which(rhs_term_names %in% possible_covariate_terms)
      if (length(remove) > 0){
        rhs_term_names <- rhs_term_names[-remove]
        parsed_rhs <- parsed_rhs[-remove]
      }
      remove <- which(rhs_term_names %in% possible_network_terms)
      if (length(remove) > 0){
        rhs_term_names <- rhs_term_names[-remove]
        parsed_rhs <- parsed_rhs[-remove]
      }

      # if there were no structural covariates in the model, then add in an
      # intercept term with a zero theta, and make sure we set our no structural
      # terms logical to TRUE
      covariate_terms_only <- FALSE
      if (length(parsed_rhs) > 0) {
        # check that the names of the statistics match those that are possible
        possible <- 1:length(rhs_term_names)
        actual <- which(rhs_term_names %in% possible_structural_terms)
        if (length(possible) != length(actual)) {
          stop(paste("the specified structural term",
                     possible_structural_terms[setdiff(possible, actual)],
                     "is not an available statistic.", sep = " "))
        }
      } else {
        covariate_terms_only <- TRUE

        # create objects by hand
        parsed_rhs <- list(edges = list(term = "edges",
             weight = 1,
             covariate = NA,
             base = NA,
             network = NA,
             threshold = 0,
             levels = NA,
             same = NA,
             method = "endogenous",
             parens_no_arg = NA,
             network_matrix_object = NA,
             num_levels = NA,
             base_index = NA))

        rhs_term_names <- "edges"

        cat("No endogenous effects were provided...\n")
      }


      # if theta is NULL, assume all ones
      if (is.null(theta) == TRUE) {
        theta <- rep(0, length(rhs_term_names))
      }
    } else {
      stop("You have specified a model with no effects, please respecify...")
    }



    # Here we are going to check and see if any of the supplied structural
    # statistics specified a $covariate field (not NA), and if they specified a
    # $base field. If they did, then we will add them to our vector.
    num_structural_stats <- 0
    statistics <- alphas <- thresholds <- thetas <- non_base_statistic <-
      theta_names <- NULL
    nodes <- list()
    if(length(parsed_rhs) > 0) {
      for (i in 1:length(parsed_rhs)) {
        term_index <- which(possible_structural_terms == rhs_term_names[i])[1]
        # if no covariate was provided then we add one statistic
        if (is.na(parsed_rhs[[i]]$covariate)) {
          num_structural_stats <- num_structural_stats + 1
          statistics <- c(statistics,term_index)
          alphas <- c(alphas, alpha[i])
          thresholds <- c(thresholds, threshold[i])
          thetas <- c(thetas, theta[i])
          nodes <- append(nodes, list(node_indicies = -1))
          non_base_statistic <- c(non_base_statistic, 0)
          theta_names <- c(theta_names, rhs_term_names[i])
        } else if (is.null(covariate_data)) {
          cat ("No covariate data was provided, but you specified a covariate for structural term:",
               rhs_term_names[i], "\n")
          # then we ignore these
          num_structural_stats <- num_structural_stats + 1
          statistics <- c(statistics,term_index)
          alphas <- c(alphas, alpha[i])
          thresholds <- c(thresholds, threshold[i])
          thetas <- c(thetas, theta[i])
          nodes <- append(nodes, list(node_indicies = -1))
          non_base_statistic <- c(non_base_statistic, 0)
          theta_names <- c(theta_names, rhs_term_names[i])
        } else {
          # check to see if the covariate name is in the column names of the
          # covariate data matrix
          ind <- which(colnames(covariate_data) == parsed_rhs[[i]]$covariate)
          if (length(ind) != 1) {
            stop (paste("invalid covariate name:",parsed_rhs[[i]]$covariate,
                        "for term:",rhs_term_names[i]))
          }

          # get the unique covariate values
          unique_covariate_values <- unique(covariate_data[,ind])

          # make sure there are atleast two categories, otherwise, this will
          # jsut be the same as the full statistic value.
          if (length(unique_covariate_values) < 2) {
            stop(paste("covariate choice of", parsed_rhs[[i]]$covariate,
                       "for term:",rhs_term_names[i],
                       "has less than two categories, please respecify!"))
          }
          # check to see if a base value was provided (so we can ignore that
          # category).
          base_value <- parsed_rhs[[i]]$base
          if (!is.na(base_value)) {
            # determine which one is the base value
            bind <- which(unique_covariate_values == base_value)[1]
            if (length(bind) != 1) {
              stop(paste("invalid base name:",base_value,
                         "for term:",rhs_term_names[i]))
            }
            unique_covariate_values <- unique_covariate_values[-bind]
          }
          # now loop over the unique covariate values, find the nodes associated
          # with them, and add everything to the list!
          for (j in 1:length(unique_covariate_values)) {
            # find the nodes associated with the unique covariate value
            node_inds <- which(covariate_data[,ind] == unique_covariate_values[j])

            if (length(node_inds) < 3) {
              stop(paste("covariate choice of", parsed_rhs[[i]]$covariate,
                         "for term:",rhs_term_names[i], "and category",
                         unique_covariate_values[j],
                         "has less than three nodes in it. All cateogries for covariates use in structural terms must have atleast three nodes. Consider using the 'base =' argument to remove this category..."))
            }

            theta_name <- paste(rhs_term_names[i], parsed_rhs[[i]]$covariate,
                                unique_covariate_values[j], sep = "_")

            num_structural_stats <- num_structural_stats + 1
            statistics <- c(statistics,term_index)
            alphas <- c(alphas, alpha[i])
            thresholds <- c(thresholds, threshold[i])
            thetas <- c(thetas, theta[i])
            nodes <- append(nodes, list(node_indicies = node_inds))
            non_base_statistic <- c(non_base_statistic, 1)
            theta_names <- c(theta_names,theta_name)
          }
        }
      }

      # check that theta and the number of statistics are equal
      if (!(length(rhs_term_names) == length(theta) |
            length(theta) == length(thetas))) {
        stop("'theta' must be the same length as the number of statistics")
      }
    }


    # for legacy implementations like Gibbs, we will need statistics in a
    # different format -- will soon be deprecated.
    stat.indx <- which(possible_structural_terms %in% rhs_term_names)
    legacy_statistics <- rep(0, num_stats)
    legacy_statistics[stat.indx] <- 1
    legacy_alphas <- rep(1, num_stats)
    legacy_thetas <- rep(0, num_stats)
    legacy_thresholds <- rep(0, num_stats)
    for (i in 1:length(rhs_term_names)) {
      legacy_alphas[which(rhs_term_names[i] == possible_structural_terms)] <- alpha[i]
      legacy_thetas[which(rhs_term_names[i] == possible_structural_terms)] <- theta[i]
      legacy_thresholds[which(rhs_term_names[i] == possible_structural_terms)] <- threshold[i]
    }
    return(list(net = net,
                statistics = statistics,
                alphas = alphas,
                thetas = thetas,
                thresholds = thresholds,
                nodes = nodes,
                non_base_statistic = non_base_statistic,
                num_structural_stats = num_structural_stats,
                legacy_statistics = legacy_statistics,
                legacy_alphas = legacy_alphas,
                legacy_thetas = legacy_thetas,
                legacy_thresholds = legacy_thresholds,
                theta_names = theta_names,
                regression_intercept = regression_intercept,
                include_intercept = include_intercept,
                covariate_terms_only = covariate_terms_only))
  }
  if(terms_to_parse == "covariate"){
    # if we are parsing covariate terms out of the formula
    remove <- which(rhs_term_names %in% possible_structural_terms)
    if (length(remove) > 0){
      rhs_term_names <- rhs_term_names[-remove]
      parsed_rhs <- parsed_rhs[-remove]
    }
    remove <- which(rhs_term_names %in% possible_network_terms)
    if (length(remove) > 0){
      rhs_term_names <- rhs_term_names[-remove]
      parsed_rhs <- parsed_rhs[-remove]
    }
    parsed_rhs <- append(parsed_rhs,list(network = net))
    return(parsed_rhs)
  }
  if(terms_to_parse == "network"){
    # if we are parsing network terms out of the formula
    remove <- which(rhs_term_names %in% possible_structural_terms)
    if (length(remove) > 0){
      rhs_term_names <- rhs_term_names[-remove]
      parsed_rhs <- parsed_rhs[-remove]
    }
    remove <- which(rhs_term_names %in% possible_covariate_terms)
    if (length(remove) > 0){
      rhs_term_names <- rhs_term_names[-remove]
      parsed_rhs <- parsed_rhs[-remove]
    }
    parsed_rhs <- append(parsed_rhs,list(network = net))
    return(parsed_rhs)
  }
}
