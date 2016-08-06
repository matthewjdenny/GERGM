calculate_additional_GOF_statistics <- function(GERGM_Object,
                                                modularity_group_memberships = NULL) {

  networks <- GERGM_Object@simulated_bounded_networks_for_GOF
  bounded_network <- GERGM_Object@bounded.network
  diag(bounded_network) <- NA
  # get the mean value of the observed network
  mean_value <- mean(bounded_network, na.rm = TRUE)

  mode <- "undirected"
  if (GERGM_Object@directed_network) {
    mode <- "directed"
  }

  if (!is.null(modularity_group_memberships)) {
    if(!is.numeric(modularity_group_memberships)) {
      stop("modularity_group_memberships must be a numeric vector indicating group memberships with minimum value 1.")
    }

    if (min(modularity_group_memberships) < 1) {
      modularity_group_memberships <- modularity_group_memberships +
        abs(min(modularity_group_memberships)) + 1
    }
  }

  mean_without_diagonal <- function(net) {
    diag(net) <- NA
    ret <- mean(net, na.rm = TRUE)
    return(ret)
  }

  get_odegrees <- function(net) {
    diag(net) <- NA
    ret <- rowSums(net, na.rm = TRUE)
    return(ret)
  }

  get_idegrees <- function(net) {
    diag(net) <- NA
    ret <- colSums(net, na.rm = TRUE)
    return(ret)
  }

  get_modularity <- function(net,
                             mode,
                             memberships) {
    net <- as.matrix(net)
    net <- igraph::graph.adjacency(net,
                                   mode = mode,
                                   weighted = TRUE,
                                   diag = FALSE)
    mod <- igraph::modularity(net, memberships, weights = igraph::E(net)$weight)
    return(mod)
  }

  
    # if we are not working in the correlation space, then intensities should
    # just be the demeaned simulated values.

    # calculate average intensity
    intensities <-  networks - mean_value

    # get the intensity of each simulated network
    simulated_intensities <- apply(intensities,3,mean_without_diagonal)

    # get the observed intensity
    observed_intensity <- mean(bounded_network - mean_value, na.rm = TRUE)

 

  # get the degree distributions on the simulated support
  simulated_odegrees <- apply(networks,3,get_odegrees)
  observed_odegrees <- get_odegrees(bounded_network)
  simulated_idegrees <- apply(networks,3,get_idegrees)
  observed_idegrees <- get_idegrees(bounded_network)

  # now calculate modularity of intensity network if group_memberships are
  # provided

  if (!is.null(modularity_group_memberships)) {

   
      # have to use raw simulated networks since weights have to be positive
      simulated_modularities <- apply(networks,
                                      3,
                                      get_modularity,
                                      mode = mode,
                                      memberships = modularity_group_memberships)
      observed_modularities <- get_modularity(bounded_network,
                                              mode,
                                              modularity_group_memberships)
    

    intensity <- c(0,observed_intensity)
    modularity <- c(0,observed_modularities)

    # now tack on stats existing to ones.
    check <- colnames(GERGM_Object@simulated_statistics_for_GOF)
    if(sum(grepl("intensity",check)) == 1 &
       sum(grepl("modularity",check) == 1)){
      ind <- which(check == "intensity")
      GERGM_Object@stats[,ind] <- intensity
      ind <- which(check == "modularity")
      GERGM_Object@stats[,ind] <- modularity
    } else if (sum(grepl("intensity",check) == 1) &
               sum(grepl("modularity",check) == 0)) {
      ind <- which(check == "intensity")
      GERGM_Object@stats[,ind] <- intensity
      GERGM_Object@stats <- cbind(GERGM_Object@stats, modularity)
    }  else {
      GERGM_Object@stats <- cbind(GERGM_Object@stats, intensity, modularity)
    }


    intensity <- simulated_intensities
    modularity <- simulated_modularities

    # need to check and see if there was already a column, if there was, then
    # put data in it, otherwise, put it in a new column.
    check <- colnames(GERGM_Object@simulated_statistics_for_GOF)
    if(sum(grepl("intensity",check)) == 1 &
       sum(grepl("modularity",check) == 1)){
      GERGM_Object@simulated_statistics_for_GOF$intensity = intensity
      GERGM_Object@simulated_statistics_for_GOF$modularity = modularity
    } else if (sum(grepl("intensity",check) == 1) &
               sum(grepl("modularity",check) == 0)) {
      GERGM_Object@simulated_statistics_for_GOF$intensity = intensity
      GERGM_Object@simulated_statistics_for_GOF <- cbind(
        GERGM_Object@simulated_statistics_for_GOF,
        modularity)
    }  else {
      GERGM_Object@simulated_statistics_for_GOF <- cbind(
        GERGM_Object@simulated_statistics_for_GOF,
        intensity,
        modularity)
    }

  } else {
    intensity <- c(0,observed_intensity)
    # now tack on stats existing to ones.
    check <- colnames(GERGM_Object@simulated_statistics_for_GOF)
    if(sum(grepl("intensity",check) == 1)){
      ind <- which(check == "intensity")
      GERGM_Object@stats[,ind] <- intensity
    }  else {
      GERGM_Object@stats <- cbind(GERGM_Object@stats, intensity)
    }

    intensity <- simulated_intensities

    check <- colnames(GERGM_Object@simulated_statistics_for_GOF)
    if(sum(grepl("intensity",check) == 1)){
      GERGM_Object@simulated_statistics_for_GOF$intensity = intensity
    }  else {
      GERGM_Object@simulated_statistics_for_GOF <- cbind(
        GERGM_Object@simulated_statistics_for_GOF,
        intensity)
    }

  }

  # add degree distribution to observed stats
  GERGM_Object@additional_stats <- list(observed_in_degrees = observed_idegrees,
                                        simulated_in_degrees = simulated_idegrees,
                                        observed_out_degrees = observed_odegrees,
                                        simulated_out_degrees = simulated_odegrees)


  return(GERGM_Object)
}
