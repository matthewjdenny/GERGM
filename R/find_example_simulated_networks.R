#' Find an example simulated network that is the minimum Frobenius  distance from the observed network.
#'
#' @param GERGM_Object The object returned by the estimation procedure using the
#' GERGM function.
#' @param observed_network The observed network used as the dependent variable
#' in the original gergm() specification.
#' @param collapse_to_correlation_space If TRUE, then all simulated networks are
#' normalized so the maximum absolute edge value is less than one. Defaults to
#' FALSE.
#' @return A simulated network (as a matrix object).
#' @export
find_example_simulated_network <- function(GERGM_Object,
                                           observed_network,
                                           collapse_to_correlation_space = F){

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
        distances[i] <- frobenius_norm(observed_network,net)
    }
    cat("\n")

    #find minimum distance
    index <- which(distances == min(distances))[1]

    cat("Minimum Frobenius Distance was:",min(distances),
        "for simulated network index number:",index, "\n")

    if (collapse_to_correlation_space) {
        net <- samples[,,index]/max(abs(samples[,,index]))
    } else {
        net <- samples[,,index]
    }

    return(net)
}
