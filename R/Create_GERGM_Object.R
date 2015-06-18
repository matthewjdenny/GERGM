# gergm object constructor (formerly gergm.object)
Create_GERGM_Object <- function(network,
                                bounded.network = network,
                                formula,
                                thetas,
                                lambda,
                                alpha,
                                together = together) {
  num.nodes <- nrow(network)
  triples <- t(combn(1:num.nodes, 3))
  h.statistics1 <- h2(network, triples = triples, statistics = rep(1, 6),
                      together = together)
  h.statistics2 <- h2(bounded.network, triples = triples,
                      statistics = rep(1, 6), together = together)
  statistics <- rbind(h.statistics1, h.statistics2)
  possible.stats <- c("out2star", "in2star", "ctriads", "recip", "ttriads",
                      "edgeweight")
  colnames(statistics) <- possible.stats
  rownames(statistics) <- c("network", "bounded.network")
  # Check whether or not the weights are NULL
  if(is.null(alpha) == TRUE){
    alpha = 0
  }
  new("gergm", network = network, bounded.network = bounded.network,
      formula = formula, stats = statistics, theta.coef = thetas,
      lambda.coef = lambda, weights = alpha)
}
