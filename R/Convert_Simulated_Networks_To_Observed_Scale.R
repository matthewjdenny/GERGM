#' Transforms simulated networks to observed scale. In general, do not use this function.
#'
#' @param GERGM_Object A GERGM object returned by the `gergm()` function. In
#' general, this function should not be used except in the case where you are
#' working with a GERGM object where the `@MCMC_output$Networks` field is still
#' on the [0,1] unconstrained space, and you wish to transform it to the
#' observed scale.
#' @return A GERGM Object
#' @export
convert_simulated_networks_to_observed_scale <- function(
  GERGM_Object){
  # determine the number of MCMC samples
  samples <- dim(GERGM_Object@MCMC_output$Networks)[3]
  num.nodes <- GERGM_Object@num_nodes
  triples = GERGM_Object@statistic_auxiliary_data$triples
  stats <- rep(1,length(GERGM_Object@stats_to_use))
  transformation_type <- GERGM_Object@transformation_type

  # figure out what kind of transformation we did, then convert back.
  printseq <- round(seq(1,samples, length.out = 11)[2:11],0)
  printcounter <- 1
  if (length(GERGM_Object@data_transformation) > 0) {
    if (GERGM_Object@is_correlation_network) {
      cat("Currently not implemented for correlation networks with covariates.")
    } else if (GERGM_Object@beta_correlation_model) {
      # if we are using a beta model for correlation networks, we can reverse
      # the transformation as follows.
      for (i in 1:samples) {
        if (i == printseq[printcounter]) {
          cat(10*printcounter,"% complete...\n", sep = "")
          printcounter <- printcounter + 1
        }
        # symmetrize incase there were any numerical imperfections
        # Symmetrize_Network(GERGM_Object@MCMC_output$Networks[,,i])
        symnet <- GERGM_Object@MCMC_output$Networks[,,i]

        # note that we stored these in the Update_Lambda_Estimates() function
        # using BZ and BZstdev as containers for mu and phi respectively.
        P <- 2*qbt(symnet, GERGM_Object@mu , GERGM_Object@phi) - 1

        # Transform P to a correlation matrix R
        R <- partials.to.correlations(P)

        # store
        GERGM_Object@MCMC_output$Networks[,,i] <- R
      }
    } else {
      for (i in 1:samples) {
        if (i == printseq[printcounter]) {
          cat(10*printcounter,"% complete...\n", sep = "")
          printcounter <- printcounter + 1
        }
        # if we did a transformation (which is the default if we are including an intercept)
        if (transformation_type == "logcauchy" |
            transformation_type == "cauchy") {
          GERGM_Object@MCMC_output$Networks[,,i] <- qst(
            GERGM_Object@MCMC_output$Networks[,,i],
            GERGM_Object@BZ,
            GERGM_Object@BZstdev,
            1)
          if (transformation_type == "logcauchy") {
            GERGM_Object@MCMC_output$Networks[,,i] <- exp(
              GERGM_Object@MCMC_output$Networks[,,i])
          }
          diag(GERGM_Object@MCMC_output$Networks[,,i]) <- 0
        }
        if (transformation_type == "lognormal" |
           transformation_type == "gaussian") {
          GERGM_Object@MCMC_output$Networks[,,i] <- qst(
            GERGM_Object@MCMC_output$Networks[,,i],
            GERGM_Object@BZ,
            GERGM_Object@BZstdev,
            Inf)
          if (transformation_type == "lognormal") {
            GERGM_Object@MCMC_output$Networks[,,i] <- exp(
              GERGM_Object@MCMC_output$Networks[,,i])
          }
          diag(GERGM_Object@MCMC_output$Networks[,,i]) <- 0
        }
        GERGM_Object@MCMC_output$Statistics[i,] <- calculate_h_statistics(
          GERGM_Object,
          GERGM_Object@statistic_auxiliary_data,
          all_weights_are_one = FALSE,
          calculate_all_statistics = TRUE,
          use_constrained_network = TRUE,
          network = GERGM_Object@MCMC_output$Networks[,,i])
      }
    }

  }else{
    if (GERGM_Object@is_correlation_network) {
      for (i in 1:samples) {
        if (i == printseq[printcounter]) {
          cat(10*printcounter,"% complete...\n", sep = "")
          printcounter <- printcounter + 1
        }
        # symmetrize incase there were any numerical imperfections
        symnet <- Symmetrize_Network(GERGM_Object@MCMC_output$Networks[,,i])
        GERGM_Object@MCMC_output$Networks[,,i] <- bounded.to.correlations(
          symnet)
      }
    } else {
      # if we did not do a transformation (only structural terms)
      cat("Currently not implemented for non-transformed networks.")
    }
  }
  return(GERGM_Object)
}
