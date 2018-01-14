#' Find an example simulated network that is the minimum Frobenius  distance from the observed network.
#'
#' @param GERGM_Object The object returned by the estimation procedure using the
#' GERGM function.
#' @param observed_network The observed network used as the dependent variable
#' in the original gergm() specification.
#' @param objective Defaults to "Frobenius", in which case the Frobenius norm
#' is used to find the simulated network that is most similar to the observed
#' network. Can also be "Likelihood", in which case the model log likelihood is
#' used to select the network to display.
#' @return A simulated network (as a matrix object).
#' @export
find_example_simulated_network <- function(GERGM_Object,
                                           observed_network,
                                           objective = c("Frobenius","Likelihood")){

  # should not need this:
  collapse_to_correlation_space = F

  possible_structural_terms <- c("out2stars",
                                 "in2stars",
                                 "ctriads",
                                 "mutual",
                                 "ttriads",
                                 "edges",
                                 "diagonal")

  # deal with objective choice
  objective <- objective[1]
  if (!(objective %in% c("Frobenius","Likelihood"))) {
    stop("objective must be one of 'Frobenius' or 'Likelihood'")
  }


    # loop over network samples and calculate similarity metric
    samples <- GERGM_Object@MCMC_output$Networks
    distances <- rep(0,dim(samples)[3])
    if (length(distances) > 50) {
        printseq <- round(seq(1,length(distances), length.out = 51)[2:51],0)
    } else {
        printseq <- 1:length(distances)
    }
    printcounter <- 1
    for (i in 1:dim(samples)[3]) {
        if (i == printseq[printcounter]) {
            cat(".")
            printcounter <- printcounter + 1
        }
        if (collapse_to_correlation_space) {
            net <- samples[,,i]/max(abs(samples[,,i]))
        } else {
            net <- samples[,,i]
        }

        if (objective == "Frobenius") {
          distances[i] <- frobenius_norm(observed_network,net)
        }
        if (objective == "Likelihood") {

          # sad <- GERGM_Object@statistic_auxiliary_data
          # hsn <- GERGM_Object@MCMC_output$Statistics[,sad$specified_statistic_indexes_in_full_statistics]


          # endogenous_contribution <- log.l(
          #   par = GERGM_Object@theta.coef[1,],
          #   alpha = GERGM_Object@weights,
          #   hsnet = hsn,
          #   ltheta = as.numeric(GERGM_Object@theta.coef[1,]),
          #   together = GERGM_Object@downweight_statistics_together,
          #   possible.stats = possible_structural_terms,
          #   GERGM_Object = GERGM_Object)

          # sub in the current simulated network
          GERGM_Object@network <- net
          if (GERGM_Object@beta_correlation_model) {
            if (net[1,1] == 0) {
              symnet <- GERGM_Object@MCMC_output$Networks[,,i]

              # note that we stored these in the Update_Lambda_Estimates() function
              # using BZ and BZstdev as containers for mu and phi respectively.
              P <- 2*qbt(symnet, GERGM_Object@mu , GERGM_Object@phi) - 1

              # Transform P to a correlation matrix R
              GERGM_Object@network <- partials.to.correlations(P)

            }
            log_likelihood <- llg.beta( par = as.numeric(GERGM_Object@lambda.coef[1,]),
                                       theta = as.numeric(GERGM_Object@theta.coef[1,]),
                                       z = GERGM_Object@data_transformation,
                                       together = GERGM_Object@downweight_statistics_together,
                                       possible.stats = possible_structural_terms,
                                 GERGM_Object = GERGM_Object)
          } else {
            log_likelihood <- llg(
              par = as.numeric(GERGM_Object@lambda.coef[1,]),
              alpha = GERGM_Object@weights,
              theta = as.numeric(GERGM_Object@theta.coef[1,]),
              z = GERGM_Object@data_transformation,
              together = GERGM_Object@downweight_statistics_together,
              possible.stats = possible_structural_terms,
              GERGM_Object = GERGM_Object)
          }



          distances[i] <- log_likelihood
        }

    }
    cat("\n")


    if (objective == "Frobenius") {
      index <- which(distances == min(distances))[1]
      cat("Minimum Frobenius Distance was:",min(distances),
          "for simulated network index number:",index, "\n")
    }
    if (objective == "Likelihood") {
      index <- which(distances == max(distances))[1]
      cat("Maximum (unnormalized) log likelihood was:",max(distances),
          "for simulated network index number:",index, "\n")
    }



    if (collapse_to_correlation_space) {
        net <- samples[,,index]/max(abs(samples[,,index]))
    } else {
        net <- samples[,,index]
    }

    return(net)
}
