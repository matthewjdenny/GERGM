# create a gergm from a formula object
Create_GERGM_Object_From_Formula <- function(object, theta.coef, together = 1, weights = NULL, 
                     transform.data = NULL, lambda.coef = NULL){
  res1 <- weighted.ergm.data(object, theta = theta.coef, 
                             alpha = weights)
  thetas <- res1$thetas
  network <- res1$net
  alphas <- res1$alphas
  statistics <- res1$statistics
  # create the network based on the transform family
  # if there are no lambda.coefficients, we assume there is no transformation
  # if there is a transformation specified, transform the observed network
  if (!is.null(lambda.coef) == TRUE){
    beta <- lambda.coef[1:(length(lambda.coef) - 1)]
    sig <- 0.01 + exp(lambda.coef[length(lambda.coef)])
    BZ <- 0
    if(is.na(dim(transform.data)[3]) == TRUE){
      BZ = BZ + beta * transform.data
    }
    if(!is.na(dim(transform.data)[3]) == TRUE){
      for (j in 1:(dim(transform.data)[3])) {
        BZ <- BZ + beta[j] * transform.data[, , j]
      }
    }
    
    bounded.network <- pst(network, BZ, sig, 1)
  }
  if (is.null(lambda.coef) == TRUE) {
    bounded.network <- network
    lambda.coef <- as.data.frame(0)
  }
  if (is.null(lambda.coef) != TRUE){
    lambda.coef <- as.data.frame(rbind(lambda.coef,NA))
    rownames(lambda.coef) <- c("est", "se")
  }
  possible.stats <- c("out2star", "in2star", "ctriads", "recip", "ttriads", 
                      "edgeweight")
  thetas <- t(as.matrix(thetas))
  thetas <- rbind(thetas, NA)
  colnames(thetas) <- possible.stats
  rownames(thetas) <- c("est", "se")
  thetas <- as.data.frame(thetas)
  
  object <- gergm.object(network = network, bounded.network = bounded.network, 
                         formula = object, thetas = thetas, 
                         lambda = lambda.coef, alpha = alphas, 
                         together = together)
}