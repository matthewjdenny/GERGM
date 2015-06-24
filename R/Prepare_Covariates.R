#Loading and pre-processing the migration data
Prepare_Covariates <- function(raw_network,
                               covariate_data,
                               covariates_to_use,
                               type_of_effect,
                               additional_network_covariates){

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

  #determine total number of covariates


  # Generate array which covaries will be transformed into
  Z <- array(0,dim=c(num_nodes,num_nodes,10))


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
