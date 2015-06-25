#Loading and pre-processing the migration data
Prepare_Network_and_Covariates <- function(raw_network,
                                    node_names = NULL,
                                    covariate_data = NULL,
                                    covariates_to_use = NULL,
                                    type_of_effect = NULL,
                                    network_covariates = NULL,
                                    network_covariate_names = NULL){

  # Calculate number of nodes in network
  num_nodes <- nrow(raw_network)

  # Check if network is square
  if (nrow(raw_network) != nrow(raw_network)) {
    stop("Network must be a square matrix or data frame!")
  }

  # If network is a data.frame, make it a matrix
  if (is.data.frame(raw_network)) {
    raw_network <- as.matrix(raw_network)
  }


  #determine whether covariates were provided and the total number of covariates
  covariates_provided <- FALSE
  # transformed covariates will have one additional slice of all ones
  num_covariates <- 1

  if(is.null(covariate_data)){
    cat("No node level covariates provided.\n")
  }else{
    covariates_provided <- TRUE
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
    covariates_provided <- TRUE
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

  # Generate array which covaries will be transformed into
  transformed_covariates <- array(0,dim=c(num_nodes,num_nodes,num_covariates))


  #1. generate sender and reciever effects
  #2. tack on any user supplied effects
  #3. standardize covariates
  #4. return list object

  #########################################################
  #Construct Z: An array of covariates which parameterize the latent space #


  Z[,,1] <- 1

  Z[,,2] <- pops07-pops06

  Z[,,3] <- popr07-popr06

  Z[,,4] <- unes07-unes06

  Z[,,5] <- uner07-uner06

  Z[,,6] <- incs07-incs06

  Z[,,7] <- incr07-incr06

  Z[,,8] <- jtemps

  Z[,,9] <- jtempr

  Z[,,10] <- statedist


  #Standardize the covariates
  for(i in 2:10){
    Z[,,i] <- (Z[,,i]-mean(c(Z[,,i])))/sd(c(Z[,,i]))
  }
  #return migr07, migr06, and Z
  return(list(Graph = migr07 - migr06, Z = Z))
}
