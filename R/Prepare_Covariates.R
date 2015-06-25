#Loading and pre-processing the migration data
Prepare_Network_and_Covariates <- function(raw_network,
                                    node_names = NULL,
                                    covariate_data = NULL,
                                    covariates_to_use = NULL,
                                    type_of_effect = NULL,
                                    network_covariates = NULL,
                                    network_covariate_names = NULL,
                                    normalization_type = c("log","division")){

  # Calculate number of nodes in network
  num_nodes <- nrow(raw_network)
  node_names <- colnames(raw_network)
  # Check if network is square
  if (nrow(raw_network) != nrow(raw_network)) {
    stop("Network must be a square matrix or data frame!")
  }

  # If network is a data.frame, make it a matrix
  if (is.data.frame(raw_network)) {
    raw_network <- as.matrix(raw_network)
  }



  #determine whether covariates were provided and the total number of covariates
  node_covariates_provided <- FALSE
  network_covariates_provided <- FALSE
  # transformed covariates will have one additional slice of all ones
  num_covariates <- 1

  if(is.null(covariate_data)){
    cat("No node level covariates provided.\n")
  }else{
    node_covariates_provided <- TRUE
    for(i in 1:length(type_of_effect)){
      if(type_of_effect[i] == "sender" | type_of_effect[i] == "reciever"){
        num_covariates <- num_covariates + 1
      }else if(type_of_effect[i] == "both"){
        num_covariates <- num_covariates + 2
      }else{
        stop("Node level covariate effect types must be one of: 'sender', 'receiver', or 'both', please respecify.")
      }
    }
    cat("You have specified",num_covariates - 1,"node level covariate effects.\n")
  }

  #determine the type and number of user provided network covariates (will not be altered)
  if(is.null(network_covariates)){
    cat("No network covariates provided.\n")
    num_additional_covars <- 0
  }else{
    network_covariates_provided <- TRUE
    if ("matrix" %in% class(network_covariates)) {
      num_covariates <- num_covariates + 1
      num_additional_covars <- 1
      cat("You have provided",num_additional_covars,"network covariates.\n")
    }
    if ("array" %in% class(network_covariates)) {
      num_covariates <- num_covariates + dim(network_covariates)[3]
      num_additional_covars <- dim(network_covariates)[3]
      cat("You have provided",num_additional_covars,"network covariates.\n")
    }
  }

  generate_covariate_effect_matrix <- function(num_nodes,
                                               node_names,
                                               covariates,
                                               covariate_column,
                                               effect_type){
    return_matrix <- matrix(0,num_nodes,num_nodes)
    if(effect_type == "sender"){
      for(j in 1:num_nodes){
        for(k in 1:num_nodes){
          if(j != k){
            row <- which(toupper(rownames(covariates)) == toupper(node_names)[j])
            return_matrix[j,k] <- covariates[row,covariate_column]
          }
        }
      }
    }
    if(effect_type == "receiver"){
      for(j in 1:num_nodes){
        for(k in 1:num_nodes){
          if(j != k){
            row <- which(toupper(rownames(covariates)) == toupper(node_names)[k])
            return_matrix[j,k] <- covariates[row,covariate_column]
          }
        }
      }
    }
    return(return_matrix)
  }

  #########################################################
  #' Construct transformed_covariates: An array of covariates which parameterize
  #' the latent space.

  # Generate array which covaries will be transformed into
  transformed_covariates <- array(0,dim=c(num_nodes,num_nodes,num_covariates))
  transformed_covariates[,,1] <- 1
  #' Set a slice counter to keep track of where we should add the covariates
  #' in the resulting covariate array.
  slice_counter <- 2

  #1. generate sender and reciever effects
  if(node_covariates_provided){

    #' Determine if type_of_effect and covariates_to_use have the same length or
    #' if no covariates_to_use is provided, that the number of columns in
    #' covariate_data is the same length.

    if(!is.null(covariates_to_use)){
      if(length(covariates_to_use) != length(type_of_effect)){
        stop("The lengths of type_of_effect and covariates_to_use are not equal, please specify vectors of equal length.")
      }
    }

    #' Loop through covariates
    for(i in 1:length(type_of_effect)){
      if(!is.null(covariates_to_use)){
        col_index <- which(tolower(colnames(covariate_data)) ==
                             tolower(covariates_to_use[i]))
      }else{
        col_index <- i
      }

      if(length(col_index) == 0){
        stop(paste("There is no matching column name in covariate_data for:",tolower(covariates_to_use[i])))
      }

      if(type_of_effect[i] == "sender"){
        add <- generate_covariate_effect_matrix(num_nodes = num_nodes,
                                                node_names = node_names,
                                                covariates = covariate_data,
                                                covariate_column = col_index,
                                                effect_type = "sender")
        transformed_covariates[,,slice_counter] <- add
        slice_counter <- slice_counter + 1
      }
      if(type_of_effect[i] == "reciever"){
        add <- generate_covariate_effect_matrix(num_nodes = num_nodes,
                                                node_names = node_names,
                                                covariates = covariate_data,
                                                covariate_column = col_index,
                                                effect_type = "receiver")
        transformed_covariates[,,slice_counter] <- add
        slice_counter <- slice_counter + 1
      }
      if(type_of_effect[i] == "both"){
        add <- generate_covariate_effect_matrix(num_nodes = num_nodes,
                                                node_names = node_names,
                                                covariates = covariate_data,
                                                covariate_column = col_index,
                                                effect_type = "sender")
        transformed_covariates[,,slice_counter] <- add
        slice_counter <- slice_counter + 1
        add <- generate_covariate_effect_matrix(num_nodes = num_nodes,
                                                node_names = node_names,
                                                covariates = covariate_data,
                                                covariate_column = col_index,
                                                effect_type = "receiver")
        transformed_covariates[,,slice_counter] <- add
        slice_counter <- slice_counter + 1
      }
    } # End of loop over node level covariates




  } # End of condition that we have node level covariates



  #2. tack on any user supplied effects
  if(network_covariates_provided){
    if(num_additional_covars == 1){
      transformed_covariates[,,slice_counter] <- network_covariates
    }else{
      for(i in 1:num_additional_covars){
        transformed_covariates[,,slice_counter] <- network_covariates[,,i]
        slice_counter <- slice_counter + 1
      }
    }
  }

  #3. standardize covariates
  for(i in 2:num_covariates){
    transformed_covariates[,,i] <- (transformed_covariates[,,i]-mean(c(transformed_covariates[,,i])))/sd(c(transformed_covariates[,,i]))
  }

  #4. return list object

  if(!node_covariates_provided & !network_covariates_provided){
    #' If no covariates were provided, then make sure the network lives on the
    #' [0,1] interval and standardize it by one of the provided methods if it
    #' does not. Then return the network and no covariates.
    if(normalization_type == "log"){
      raw_network <- raw_network
      network <-
    }

    return(network)
  }else{
    #' If covariates were provided, put the network and the covariates in a list
    #' object and return it.


    return(list(network = raw_network, transformed_covariates = transformed_covariates))
  }

} # End of function definition.
