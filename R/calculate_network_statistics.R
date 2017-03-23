#' @title A Function to estimate calculate (weighted) network statistics
#' @description Calculates out2stars, in2stars, ctriads, mutual, ttriads, and
#' edges statistics (with or without exponential down weighting) for a
#' real-valued network.
#'
#' @param network A square numeric matrix (sociomatrix or adjacency matrix)
#' representing the network.
#' @param weights If you wish to provide your own weights, you must provide a
#' vector of length 6 with terms corresponding to out2stars, in2stars, ctriads,
#' mutual, ttriads, edges in that order.
#' @param downweight_statistics_together Logical indicating whether exponential
#' down weighting should be done together or separately. Defaults to TRUE.
#' @return A gergm object containing parameter estimates.
#' @param include_diagonal Logical indicating whether the diagonal should be
#' included in the network statistics. Defaults to FALSE.
#' @export
calculate_network_statistics <- function(
  network,
  weights = NULL,
  downweight_statistics_together = TRUE,
  include_diagonal = FALSE) {


  # hard coded possible stats
  possible_stats <- c("out2stars",
                      "in2stars",
                      "ctriads",
                      "mutual",
                      "ttriads",
                      "edges",
                      "diagonal")

  #if no weights are provided, make them all 1's
  if (is.null(weights)){
    weights <- rep(1, length(possible_stats))
  } else if (length(weights) != length(possible_stats)) {
    stop(paste("If you wish to provide your own weights, you must provide a vector of length",
               length(possible_stats),"with terms corresponding to",
               paste0(possible_stats,collapse = ", "), "in that order."))
  }

  # makre sure it is square
  if (nrow(network) != ncol(network)) {
    stop("The argument 'network' must be a square numeric matrix.")
  }

  # make sure the diagonal is zeros
  diag(network) <- 0
  # get the triples to pass in
  num.nodes <- nrow(network)
  if (include_diagonal) {
    triples <- t(combn(1:num.nodes, 3))
    # now add in all triples that include the diagonal as pairs will just be
    # captured in the "diagonal" term. These will just be of the form (i,i,j)
    # there will therefore be
    num_to_add <- num.nodes * (num.nodes - 1)
    add <- matrix(0,nrow = num_to_add,ncol = 3)
    add_counter <- 1
    for (i in 1:num.nodes) {
      for (j in 1:num.nodes) {
        if (i != j) {
          add[add_counter,] <- c(i,i,j)
          add_counter <- add_counter + 1
        }
      }
    }
    triples <- rbind(triples, add)
  } else {
    triples <- t(combn(1:num.nodes, 3))
  }

  # calculate the statistics
  statistics <- h2(network,
                   triples = triples,
                   statistics = rep(1, length(possible_stats)),
                   alphas = weights,
                   together = downweight_statistics_together)

  #put everything in a data.frame
  to_return <- data.frame(out2stars = statistics[1],
                          in2stars  = statistics[2],
                          ctriads = statistics[3],
                          mutual = statistics[4],
                          ttriads = statistics[5],
                          edges = statistics[6],
                          diagonal = statistics[7],
                          stringsAsFactors = FALSE)

  return(to_return)
}
