# gergm object constructor (formerly gergm.object)
Create_GERGM_Object <- function(network,
                                bounded.network = network,
                                formula,
                                thetas,
                                lambda,
                                alpha,
                                together = together,
                                possible.stats,
                                thresholds = NULL) {
  if(is.null(thresholds)){
    thresholds <- rep(0, length(possible.stats))
  }
  num.nodes <- nrow(network)
  # Check whether or not the weights are NULL
  if(is.null(alpha) == TRUE){
    alpha = 0
  }
  new("gergm", network = network,
      bounded.network = bounded.network,
      formula = formula,
      theta.coef = thetas,
      lambda.coef = lambda,
      weights = alpha,
      num_nodes = num.nodes,
      thresholds = thresholds)
}
