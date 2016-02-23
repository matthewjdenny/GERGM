##Transformation functions for correlation matrices
##James D. Wilson 2/18/15
#-----------------------------------------------------
correlations.to.partials <- function(correlations){
  #INPUT:
  #      correlations - d x d matrix whose (i,j)th entry is the
  #                     correlation between i and j
  #                     NOTE: diagonal should be all set to 1
  #OUTPUT:
  #       partials - d x d matrix whose (i,j)th entry is the partial
  #                  correlation between i and j given all variables
  #                  in between i and j

  #Initialize
  d <- dim(correlations)[1]
  partials <- matrix(1, d, d)

  #NOTE: partials should be symmetric. If not, send an error
  if (isSymmetric(correlations) == FALSE) {
    stop("The correlation matrix must be symmetric.")
  }

  #if(sum(eigen(correlations)$values < 0) > 0){
  #  stop("The correlation matrix must be positive definite.")
  #}

  #Set the super- and sub- diagonal of the partials matrix
  #super-diagonal
  diag(partials[-nrow(partials), -1]) <- diag(correlations[-nrow(correlations), -1])

  #sub-diagonal
  diag(partials[-1, -ncol(partials)]) <- diag(correlations[-1, -ncol(correlations)])

  #Recursively obtain the remaining correlations using the transformation
  #given in the Harry Joe paper.

  for (k in 2:(d - 1)) {
    for (i in 1:(d - k)) {

      R2 <- correlations[(i + 1):(i + k - 1), (i + 1):(i + k - 1)]
      r1 <- correlations[i, (i + 1):(i + k - 1)]
      r3 <- correlations[i + k, (i + 1):(i + k - 1)]

      D <- sqrt((1 - t(r1) %*% solve(R2) %*% r1) *
                  (1 - t(r3) %*% solve(R2) %*% r3))

      partials[i, i + k] <- as.numeric((correlations[i, i + k] -
                                          t(r1) %*% solve(R2) %*% r3) / D)

      partials[i + k, i] <- partials[i, i + k]
    }
  }
  return(partials)
}
#-----------------------------------------------------
partials.to.correlations <- function(partials){
  #INPUT:
  #       partials - d x d matrix whose (i,j)th entry is the partial
  #                  correlation between i and j given all variables
  #                  in between i and j
  #OUTPUT:
  #       correlations - d x d matrix whose (i,j)th entry is the
  #                     correlation between i and j
  #                     NOTE: diagonal should be all set to 1


  #Initialize
  d <- dim(partials)[1]
  correlations <- matrix(1, d, d)

  #NOTE: partials should be symmetric. If not, send an error
  if (isSymmetric(partials) == FALSE) {
    stop("The partial correlation matrix must be symmetric.")
  }

  #super-diagonal
  diag(correlations[-nrow(correlations), -1]) <- diag(partials[-nrow(partials), -1])

  #sub-diagonal
  diag(correlations[-1, -ncol(correlations)]) <- diag(partials[-1, -ncol(partials)])

  #NOTE: the above should be symmetric. If not, send an error
  if (isSymmetric(correlations) == FALSE) {
    stop("The correlation matrix must be symmetric. Check the super- and
         sub- diagonals of object partials.")
  }

  #Recursively obtain the remaining correlations using the transformation
  #given in the Harry Joe paper.

  for(k in 2 : (d - 1)){
    for(i in 1 : (d - k)){

      R2 <- correlations[(i + 1): (i + k - 1), (i + 1): (i + k - 1)]
      r1 <- correlations[i, (i + 1) : (i + k - 1)]
      r3 <- correlations[i + k, (i + 1): (i + k - 1)]

      D <- sqrt((1 - t(r1) %*% solve(R2) %*% r1)* (1 - t(r3) %*% solve(R2) %*% r3))
      correlations[i, i + k] <- t(r1) %*% solve(R2) %*% r3 + partials[i, i + k]*D
      correlations[i + k, i] <- correlations[i, i + k]
    }
  }
  return(correlations)
  }

#-----------------------------------------------------
#Function to transform a correlation network to a bounded network
#via the beta transform
pbt <- function(correlations, mu, phi){
  #Note: mu must be n choose 2 in length, and phi a single parameter
  #Takes a symmetric matrix correlations and transforms this to an undirected
  #network X on [0,1] where X is given as an edge list.
  d <- dim(correlations)[1] #number of nodes
  X <- diag(d)
  #transform to partial correlation network first
  P <- correlations.to.partials(correlations)
  shape1 <- mu * phi
  shape2 <- (1 - mu)*phi
  #transform to bounded network now
  probs <- pbeta(P[lower.tri(P, diag = FALSE)],
                 shape1 = shape1, shape2 = shape2)
  X[lower.tri(X)] <- probs
  #make symmetric
  X <- X + t(X)
  diag(X) <- 0
  return(X)
}
#-----------------------------------------------------
#Function to calculate Beta density function evaluated at the observed network
dbt <- function(x, mu, phi){
  #Note: x and mu must be n choose 2 in length, and phi a single parameter
  #calculate parameters first
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  density <- dbeta(x, shape1 = shape1, shape2 = shape2)
  return(density)
}
#-----------------------------------------------------
#Function to calculate the inverse of a Beta distribution function
#evaluated at the observed (constrained) network X on [0,1]
qbt <- function(X, mu, phi){
  #Note: mu must be an n choose 2 vector, where n is the number of nodes in X

  n <- dim(X)[1] #number of nodes
  P <- diag(n)

  shape1 <- mu*phi
  shape2 <- (1 - mu) * phi
  quants <- qbeta(X[lower.tri(X, diag = FALSE)],
                  shape1 = shape1, shape2 = shape2)
  P[lower.tri(P)] <- quants
  #make symmetric
  P <- P + t(P)
  diag(P) <- 1
  return(P)
}
#-----------------------------------------------------
#Logistic function: to get mu from regression coefficients beta
logistic <- function(beta, z){
  #INPUT:
  # beta - vector of regression coefficients of length k
  # z - k-dimensional array of n x n edge valued statistics (must be an array! and first dimension
  # should be all 1s for an intercept)
  #OUTPUT:
  # mu - a n choose 2 vector of edge means
  BZ <- 0

  for (i in 1 : (dim(z)[3])) {
    BZ <- BZ + beta[i]  * z[, , i]
  }

  BZ <- BZ[lower.tri(BZ)]
  mu <- exp(BZ) / (1 + exp(BZ))
  return(mu)
}

#-----------------------------------------------------
#Function for optimizing mean (mu) and phi for the Beta transform
llg.beta <- function(par,
                     theta,
                     z,
                     together = together,
                     possible.stats = possible.stats,
                     GERGM_Object = GERGM_Object) {
  # log likelihood for unbounded network with g function
  statistics <- GERGM_Object@stats_to_use
  alphas <- GERGM_Object@weights
  net <- GERGM_Object@network #call original correlation network

  #get parameters
  beta <- par[1:(length(par) - 1)] #beta should be k + 1 in length (including at least an intercept)
  phi <- par[length(par)] #alpha should be 1 in length

  #From beta, get mean edge weights
  #This is supposing that z is an array of n x n matrices corresponding to each statistic
  #This also assumes that z is at least of dimension 1 (to include an intercept)

  mu <- logistic(beta, z)

  #create bounded network via beta transform
  net2 <- pbt(net, mu, phi)

  num.nodes <- nrow(net2)
  triples <- t(combn(1:num.nodes, 3))

  #calculate log likelihood over which we are optimizing
  log.li <- rbind(theta) %*%
    h2(net2,
       triples = triples,
       statistics = statistics,
       alphas = alphas,
       together = together) +
    2*sum(log(dbt((net[upper.tri(net)]+1)/2, mu, phi)))

  return(as.numeric(log.li))
}
