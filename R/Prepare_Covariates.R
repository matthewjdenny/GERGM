#Loading and pre-processing the migration data
Prepare_Covariates <- function(raw_network,
                               covariate_data,
                               type_of_effect){

  #1. generate sender and reciever effects
  #2. tack on any user supplied effects
  #3. standardize covariates
  #4. return list object

  #########################################################
  #Construct Z: An array of covariates which parameterize the latent space #

  Z <- array(0,dim=c(nrow(migr05),nrow(migr05),10))
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
