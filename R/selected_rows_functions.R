generate_selected_rows <- function(pairs,
                                   triples,
                                   endogenous_statistic_node_sets,
                                   full_statistics,
                                   non_base_indicator) {
  # indicate
  # pairs_statistics <- c(4,6)
  triples_statistics <- c(1,2,3,5)

  # find the index of the first non-base statistic

  # determine whether we are evn making node sets. If not, then just return
  # stuff and move on
  creating_node_sets <- FALSE
  if (sum(non_base_indicator > 0)) {
    creating_node_sets <- TRUE
    # find the index of the first non-base statistic
    first_non_base <- min(which(non_base_indicator == 1)) - 1
  }

  if (!creating_node_sets) {
    # set these to two so that we avoid uint conversion issues.
    save_statistics_selected_rows_matrix <- matrix(2,
                                                   nrow = 2,
                                                   ncol = length(full_statistics))
    rows_to_use <- rep(2, length(full_statistics))
  } else {
    # determine which entries in the
    max_length <- 0
    full_stat_counter <- 1
    indicies <- vector(mode = "list",
                       length = length(endogenous_statistic_node_sets))
    for (i in 1:length(endogenous_statistic_node_sets)) {
      node_indicies <- endogenous_statistic_node_sets[[i]]
      # if we are dealing with actual nodes
      if (length(node_indicies) > 2) {
        if (full_statistics[first_non_base + full_stat_counter] %in%
            triples_statistics) {
          indicies[[i]] <- get_triples_rows(triples, node_indicies)
        } else {
          indicies[[i]] <- get_pairs_rows(pairs, node_indicies)
        }
        max_length <- max(max_length,length(indicies[[i]]))
        full_stat_counter <- full_stat_counter + 1
      }
    }

    # now we generate the rows matrix
    save_statistics_selected_rows_matrix <- matrix(0, nrow = max(2,max_length),
                                                   ncol = length(full_statistics))
    rows_to_use <- rep(0, length(full_statistics))

    # remove zero length indicies entries
    to_rm <- which(sapply(indicies, length) == 0)
    indicies <- indicies[-to_rm]

    # populate the matrix
    counter <- 1
    for (i in 1:length(full_statistics)) {
      if(non_base_indicator[i] > 0) {
        l <- length(indicies[[counter]])
        save_statistics_selected_rows_matrix[1:l,i] <- indicies[[counter]]
        rows_to_use[i] <- l
        counter <- counter + 1
      }
    }
  }

  return(list(save_statistics_selected_rows_matrix = save_statistics_selected_rows_matrix,
              rows_to_use = rows_to_use))
}


get_triples_rows <- function(triples,
                           node_indicies) {
  rows <- which(triples[,1] %in% node_indicies &
                triples[,2] %in% node_indicies &
                triples[,3] %in% node_indicies)
  rows <- rows[order(rows, decreasing = FALSE)]
  return(rows)
}

get_pairs_rows <- function(pairs,
                            node_indicies) {
  rows <- which(pairs[,1] %in% node_indicies &
                pairs[,2] %in% node_indicies)
  rows <- rows[order(rows, decreasing = FALSE)]
  return(rows)
}

# stochastic_MH_proportion <- GERGM_Object@stochastic_MH_proportion
generate_stochastic_MH_triples_pairs <- function(stochastic_MH_proportion,
                                                 use_stochastic_MH,
                                                 triples,
                                                 pairs,
                                                 samples = 100) {

  if (use_stochastic_MH) {
    # first deal with triples
    ntrip <- nrow(triples)
    # make sure it is always atleast two so that we get back a matrix
    ntrip <- max(ceiling(ntrip * stochastic_MH_proportion), 2)
    random_triad_samples <- matrix(0, nrow = ntrip, ncol = samples)
    for (i in 1:samples) {
      cur <- sample(x = 1:nrow(triples), size = ntrip, replace = FALSE)
      cur <- cur[order(cur,decreasing = FALSE)]
      random_triad_samples[,i] <- cur
    }

    # now deal with pairs
    np <- nrow(pairs)
    # make sure it is always atleast two so that we get back a matrix
    np <- max(ceiling(np * stochastic_MH_proportion), 2)
    random_dyad_samples <- matrix(0, nrow = np, ncol = samples)
    for (i in 1:samples) {
      cur <- sample(x = 1:nrow(pairs), size = np, replace = FALSE)
      cur <- cur[order(cur,decreasing = FALSE)]
      random_dyad_samples[,i] <- cur
    }
  } else {
    # if we are not using stochastic MH, then populate a matrix of all zeros
    random_triad_samples <- matrix(1, nrow = 2, ncol = samples)
    random_dyad_samples <- matrix(1, nrow = 2, ncol = samples)
  }

  # now generate a list of subsets of the pairs and triples matrices
  random_triples <- vector(mode = "list", length = samples)
  random_pairs <- vector(mode = "list", length = samples)

  for (i in 1:samples) {
    random_triples[[i]] <- triples[random_triad_samples[,i],] - 1
    random_pairs[[i]] <- pairs[random_dyad_samples[,i],] - 1
  }

  return(list(random_triad_samples = random_triad_samples,
              random_dyad_samples = random_dyad_samples,
              random_triples = random_triples,
              random_pairs = random_pairs))
}


prepare_statistic_auxiliary_data <- function(GERGM_Object) {

  num_nodes <- GERGM_Object@num_nodes
  triples <- t(combn(1:num_nodes, 3))
  pairs <- t(combn(1:num_nodes, 2))
  endogenous_statistic_node_sets <- GERGM_Object@endogenous_statistic_node_sets

  # determine which statistics we are using when calculating
  non_base_stats <- GERGM_Object@stats_to_use[
    which(GERGM_Object@non_base_statistic_indicator == 1)]
  if (length(non_base_stats) > 0) {
    if (GERGM_Object@directed_network) {
      full_statistics <- c(1:6,non_base_stats)
      non_base_indicator <- c(rep(0,6), rep(1,length(non_base_stats)))
      base_statistics_to_save <- 1:6

    } else {
      full_statistics <- c(2,5,6,non_base_stats)
      non_base_indicator <- c(rep(0,3), rep(1,length(non_base_stats)))
      base_statistics_to_save <- c(2,5,6)
    }
  } else {
    if (GERGM_Object@directed_network) {
      full_statistics <- c(1:6)
      non_base_indicator <- rep(0,6)
      base_statistics_to_save <- 1:6
    } else {
      full_statistics <- c(2,5,6)
      non_base_indicator <- rep(0,3)
      base_statistics_to_save <- c(2,5,6)
    }
  }

  base_statistic_alphas <- rep(1,length(base_statistics_to_save))
  print(base_statistic_alphas)
  # need to set the base statistic alphas to less than one if they were specified
  # to have a value less than one.
  bs <- GERGM_Object@stats_to_use[
    which(GERGM_Object@non_base_statistic_indicator == 0)]
  if (length(bs) > 0) {
    for (i in 1:length(bs)) {
      ind <- which(base_statistics_to_save == bs[i])
      base_statistic_alphas[ind] <- GERGM_Object@weights[i]
    }
  }
  # generate rows to use
  result <- generate_selected_rows(pairs,
                                   triples,
                                   endogenous_statistic_node_sets,
                                   full_statistics,
                                   non_base_indicator)

  selected_rows_matrix <- result$save_statistics_selected_rows_matrix
  rows_to_use <- result$rows_to_use


  # now deal with the case where we simply want to work with the actual statistics
  # and not add in anything.
  base_statistics_to_save2 <- GERGM_Object@stats_to_use[which(
    GERGM_Object@non_base_statistic_indicator == 0
  )]
  base_statistic_alphas2 <- GERGM_Object@weights

  result <- generate_selected_rows(pairs,
                                   triples,
                                   endogenous_statistic_node_sets,
                                   GERGM_Object@stats_to_use,
                                   GERGM_Object@non_base_statistic_indicator)

  selected_rows_matrix2 <- result$save_statistics_selected_rows_matrix
  rows_to_use2 <- result$rows_to_use


  # figure out which entries in the full stats correspond to the stats that
  # the user specified
  ind1 <- GERGM_Object@stats_to_use[
    which(GERGM_Object@non_base_statistic_indicator == 0)]

  # need to map directed to undirected stats
  if (!GERGM_Object@directed_network) {
    stat_mapping <- GERGM_Object@possible_endogenous_statistic_indices
    for (l in 1:length(ind1)) {
      ind1[l] <- which(stat_mapping == ind1[l])[1]
    }
  }

  ind2 <- which(GERGM_Object@non_base_statistic_indicator == 1) +
    length(base_statistics_to_save) -
    length(ind1)

  specified_statistic_indexes_in_full_statistics <- c(ind1,ind2)

  return(list(full_statistics = full_statistics,
              full_non_base_indicator = non_base_indicator,
              full_base_statistics_to_save = base_statistics_to_save,
              full_base_statistic_alphas = base_statistic_alphas,
              full_selected_rows_matrix = selected_rows_matrix,
              full_rows_to_use = rows_to_use,
              specified_statistic_indexes_in_full_statistics =
                specified_statistic_indexes_in_full_statistics,
              specified_base_statistics_to_save = base_statistics_to_save2,
              specified_base_statistic_alphas = base_statistic_alphas2,
              specified_selected_rows_matrix = selected_rows_matrix2,
              specified_rows_to_use = rows_to_use2))
}

