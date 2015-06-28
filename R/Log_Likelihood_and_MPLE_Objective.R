# log likelihood
log.l <- function(thetas,
                  alpha,
                  hsnet,
                  ltheta,
                  together = together,
                  possible.stats,
                  GERGM_Object = GERGM_Object) {
  #turn dataframe into matrix
  hsnet <- as.matrix(hsnet)
  if (nrow(hsnet) == 1) {
    theta <- thetas[1:nrow(hsnet)]
    z <- hsnet * (theta - ltheta)
  }
  if (nrow(hsnet) > 1) {
    theta <- thetas[1:ncol(hsnet)]
    #print(str(hsnet))
    z <- hsnet %*% (theta - ltheta)
  }

  temp <- h(possible.stats,
            alpha,
            theta = theta,
            together = together,
            GERGM_Object)[1, ]
  return(rbind(theta) %*% temp - max(z) - log(sum(exp(z - max(z)))))
}

llg <- function(par,
                alpha,
                theta,
                z,
                together = together,
                possible.stats = possible.stats,
                GERGM_Object = GERGM_Object) {
  # log likelihood for unbounded network with g function
  statistics <- GERGM_Object@stats_to_use
  alphas <- GERGM_Object@weights
  net <- GERGM_Object@network
  beta <- par[1:(length(par) - 1)]
  sig <- 0.01 + exp(par[length(par)])
  BZ <- 0
  for (i in 1:(dim(z)[3])) {
    BZ <- BZ + beta[i] * z[, , i]
  }
  net2 <- pst(net, BZ, sig, 1)
  num.nodes <- nrow(net2)
  triples <- t(combn(1:num.nodes, 3))
  temp <- h2(net2,
             triples = triples,
             statistics = statistics,
             alphas = alphas,
             together = together)
  log.li <- rbind(theta) %*%
    h2(net2,
       triples = triples,
       statistics = statistics,
       alphas = alphas,
       together = together) +
    sum(log(dst(net[upper.tri(net)], BZ[upper.tri(net)], sig, 1))) +
    sum(log(dst(net[lower.tri(net)], BZ[lower.tri(net)], sig, 1)))
  return(as.numeric(log.li))
}

# maximum pseudo-likelihood estimates
mple <- function(net, statistics, directed) {
  xy <- net2xy(net, statistics, directed = directed)
  x <- xy$x
  y <- xy$y
  est <- coef(lm(y ~ x - 1))
  ests <- optim(par = est, pl, y = y, x = x, method = "BFGS",
                hessian = TRUE,control = list(fnscale = -1, trace = 6))
  return(ests)
}
