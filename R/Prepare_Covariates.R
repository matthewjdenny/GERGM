#' A Function to prepare network and covariate data for use in the gergm() function.
#'
#' @param raw_network A square matrix or data frame object containing the network we wish to estimate the gergm parameters for. If not providing node_names but using covariate_data, it should have identical row and column names that correspond to the rownames in the covariate_data matrix so as to facilitate matching of the node level covariates to rows/columns in the adjacency matrix.
#' @param node_names An optional character vector of node names that will be used to match the row observations of covariate_data to the corresponding nodes in the network. If specified, these will override existing row and column names in the raw_network object.
#' @param covariate_data A matrix or data frame containing node level covariates the user wished to transform into sender or reciever effects. It must have row names that match every entry in node_names or colnames(raw_network) if node_names = NULL, should have descriptive column names. More rows may be provided than will be used in the analysis and additional columns may also be provided if specifying covariates_to_use. If left NULL, then no sender or reciever effects will be added.
#' @param covariates_to_use A vector of column names that correspond to a subset of the columns (and their exact column names) in the covariate_data. This is an easy way to only specify a subset of covariates in a node level dataset. If covariates_to_use = NULL, then all columns will be used.
#' @param type_of_effect A character vector specifying the type of effec(s) to be included in the resulting transformed covariates array for each covariate. Each entry may be one of "sender", "receiver", or "both", correponding to the relevant effects being included for each covariate. If left NULL, then no covariate effects will be specified.
#' @param network_covariates A matrix (or array, for multiple effects) containing additional network covariates the user wishes to include in the model that are not sender or receiver effects (such as geographic distance, or group co-membership). The user may specify arbitrary effects in this matrix (array) to be included in the estimation.
#' @param network_covariate_names A vector of names for each of the slices (or the single matrix) of network_covariates. If no names are specified, the each of these covariates will be automatically named network_covariate_1, network_covariate_2, ...
#' @param normalization_type If only a raw_network is provided then, the function will automatically check to determine if all edges fall in the [0,1] interval. If edges are determined to fall outside of this interval, then a trasformation onto the interval may be specified. If "division" is selected, then the data will have a value added to them such that the minimum value is atleast zero (if necessary) and then all edge values will be divided by the maximum to ensure that the maximum value is in [0,1]. If "log" is selected, then the data will have a value added to them such that the minimum value is atleast zero (if necessary), then 1 will be added to all edge values before they are logged and then divided by the largest value, again ensuring that the resulting network is on [0,1]. Defaults to "log" and need not be set to NULL if providing covariates as it will be ignored.
#' @return If only the raw_network is provided then, this function returns the normalized bounded network (if normalization is required). If covariates are provided, then it returns a list object with a $network field that can be assigned to the network object to be included in the formula for the gergm() function for estimation and a $transformed_covariates field that can be assigned to the data_transformation parameter of the gergm() function.
#' @export
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
  if(is.null(node_names)){
    node_names <- colnames(raw_network)
  }

  # Make sure that the diagnoal terms in the network are zeros
  diag(raw_network) <- 0

  # Check if network is square
  if (nrow(raw_network) != nrow(raw_network)) {
    stop("Network must be a square matrix or data frame!")
  }

  # Just to make sure the code does not break!
  if(is.null(normalization_type)){
    normalization_type <- "log"
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




  #3. return list object

  if(!node_covariates_provided & !network_covariates_provided){
    #' If no covariates were provided, then make sure the network lives on the
    #' [0,1] interval and standardize it by one of the provided methods if it
    #' does not. Then return the network and no covariates.

    if(min(raw_network) < 0){
      raw_network <- raw_network - min(raw_network)
    }
    if(max(raw_network) > 1){
      if(normalization_type[1] == "log"){
        raw_network <- raw_network + 1
        network <- log(raw_network)
        network <- network/max(network)
      }

      if(normalization_type[1] == "division"){
        network <- raw_network/max(raw_network)
      }

      diag(network) <- 0
    }else{
      network <- raw_network
    }

    return(network)
  }else{
    #4. standardize covariates
    for(i in 2:num_covariates){
      transformed_covariates[,,i] <- (transformed_covariates[,,i]-mean(c(transformed_covariates[,,i])))/sd(c(transformed_covariates[,,i]))
    }

    #' If covariates were provided, put the network and the covariates in a list
    #' object and return it. First we want to build up a list of covariate names
    #' which we will use as the slice names for the transformed covariate object
    slice_names <- "intercept"
    slice_counter <- 2
    if(node_covariates_provided){
      for(i in 1:length(type_of_effect)){
        if(type_of_effect[i] == "sender"){
          if(!is.null(covariates_to_use)){
            slice_names <- c(slice_names, paste(covariates_to_use[i],"_sender",sep = ""))
          }else{
            slice_names <- c(slice_names,paste(colnames(covariate_data)[i],"_sender",sep = ""))
          }
        }
        if(type_of_effect[i] == "reciever"){
          if(!is.null(covariates_to_use)){
            slice_names <- c(slice_names,paste(covariates_to_use[i],"_receiver",sep = ""))
          }else{
            slice_names <- c(slice_names,paste(colnames(covariate_data)[i],"_reciever",sep = ""))
          }

        }
        if(type_of_effect[i] == "both"){
          if(!is.null(covariates_to_use)){
            slice_names <- c(slice_names,paste(covariates_to_use[i],"_sender",sep = ""))
          }else{
            slice_names <- c(slice_names,paste(colnames(covariate_data)[i],"_sender",sep = ""))
          }
          if(!is.null(covariates_to_use)){
            slice_names <- c(slice_names,paste(covariates_to_use[i],"_receiver",sep = ""))
          }else{
            slice_names <- c(slice_names,paste(colnames(covariate_data)[i],"_reciever",sep = ""))
          }
        }
      }
    }#node covariates provided conditional

    if(network_covariates_provided){
      if(!is.null(network_covariate_names)){
        slice_names <- c(slice_names,network_covariate_names)
      }else{
        temp <- NULL
        for(i in 1:num_additional_covars){
          temp <- c(temp,paste("network_covariate_",i,sep = ""))
        }
        slice_names <- c(slice_names,temp)
      }
    }

    #assign the dimnames to the array object
    dimnames(transformed_covariates) <- list(node_names,node_names,slice_names)

    return(list(network = raw_network, transformed_covariates = transformed_covariates))
  }

} # End of function definition.
