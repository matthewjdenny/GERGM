Convert_Simulated_Networks_To_Observed_Scale <- function(
  GERGM_Object,
  transformation_type){

  # determine hte number of MCMC samples
  samples <- dim(GERGM_Object@MCMC_Output$Networks)[3]
  num.nodes <- GERGM_Object@num_nodes
  triples = t(combn(1:num.nodes, 3))
  stats <- rep(1,length(GERGM_Object@stats_to_use))

  # figure out what kind of transformation we did, then convert back.
  if (length(GERGM_Object@data_transformation) > 0) {
    for(i in 1:samples){
      # if we did a transformation (which is the default if we are including an intercept)
      if(transformation_type == "logcauchy" | transformation_type == "cauchy"){
        GERGM_Object@MCMC_Output$Networks[,,i] <- qst(
          GERGM_Object@MCMC_Output$Networks[,,i],
          GERGM_Object@BZ,
          GERGM_Object@BZstdev,
          1)
      }
      if(transformation_type == "lognormal" | transformation_type == "gaussian"){
        GERGM_Object@MCMC_Output$Networks[,,i] <- qst(
          GERGM_Object@MCMC_Output$Networks[,,i],
          GERGM_Object@BZ,
          GERGM_Object@BZstdev,
          Inf)
      }
      GERGM_Object@MCMC_Output$Statistics[i,] <- h2(
        net = GERGM_Object@MCMC_Output$Networks[,,i],
        triples = triples,
        statistics = stats,
        alphas = NULL,
        together = 1,
        directed = TRUE)
    }
  }else{
    # if we did not do a transformation (only structural terms)
    cat("Currently not implemented for non-transformed networks.")
  }
  return(GERGM_Object)
}
