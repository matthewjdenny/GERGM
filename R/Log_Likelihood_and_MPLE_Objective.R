# log likelihood
log.l <- function(thetas, formula, alpha, hsnet, ltheta, together = together, possible.stats) {
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

  temp <- h(formula, possible.stats, alpha, theta = theta, together = together)[1, ]
  return(rbind(theta) %*% temp - max(z) - log(sum(exp(z - max(z)))))

  #theta <- par[1:ncol(hsnet)]
  #z <- hsnet%*%(theta-ltheta)
  #rbind(theta)%*%h(net,triples)-max(z)-log(sum(exp(z-max(z))))
  #return(rbind(theta) %*% temp - max(z) - log(mean(exp(z - max(z)))))
}

llg <- function(par, formula, alpha, theta, z, together = together, possible.stats = possible.stats) {
  # log likelihood for unbounded network with g function
  res1 <- Parse_Formula_Object(formula, possible.stats, theta = theta, alpha = alpha)
  statistics <- res1$statistics
  alphas <- res1$alphas
  net <- res1$net
  beta <- par[1:(length(par) - 1)]
  sig <- 0.01 + exp(par[length(par)])
  BZ <- 0
  #print(gpar)
  for (i in 1:(dim(z)[3])) {
    BZ <- BZ + beta[i] * z[, , i]
  }
  #print(beta) #added this for check
  net2 <- pst(net, BZ, sig, 1)
  num.nodes <- nrow(net2)
  #print(net.new)
  triples <- t(combn(1:num.nodes, 3))
  log.li <- rbind(theta) %*%
    h2(net2, triples = triples, statistics = statistics, alphas = alphas,
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
