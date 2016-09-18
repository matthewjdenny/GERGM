
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

  if(GERGM_Object@is_correlation_network){
    #this will calculate the h statistics on the original network as desired
    temp <- h.corr(possible.stats,
                   alpha,
                   theta = theta,
                   together = together,
                   GERGM_Object)[1, ]
  }else{
    temp <- calculate_h_statistics(
      GERGM_Object,
      GERGM_Object@statistic_auxiliary_data,
      all_weights_are_one = FALSE,
      calculate_all_statistics = FALSE,
      use_constrained_network = TRUE)
  }
  ret <- rbind(theta) %*% temp - max(z) - log(sum(exp(z - max(z))))
  return(ret)
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
  transformation_type <- GERGM_Object@transformation_type
  if(transformation_type == "logcauchy"){
    net2 <- pst(log(net), BZ, sig, 1)
    last_term <- sum(log(dst(log(net[upper.tri(net)]), BZ[upper.tri(net)], sig, 1))) +
      sum(log(dst(log(net[lower.tri(net)]), BZ[lower.tri(net)], sig, 1)))
  }
  if( transformation_type == "cauchy"){
    net2 <- pst(net, BZ, sig, 1)
    last_term <- sum(log(dst(net[upper.tri(net)], BZ[upper.tri(net)], sig, 1))) +
      sum(log(dst(net[lower.tri(net)], BZ[lower.tri(net)], sig, 1)))
  }
  if(transformation_type == "lognormal"){
    net2 <- pst(log(net), BZ, sig, Inf)
    last_term <- sum(log(dst(log(net[upper.tri(net)]), BZ[upper.tri(net)], sig, Inf))) +
      sum(log(dst(log(net[lower.tri(net)]), BZ[lower.tri(net)], sig, Inf)))
  }
  if( transformation_type == "gaussian"){
    net2 <- pst(net, BZ, sig, Inf)
    last_term <- sum(log(dst(net[upper.tri(net)], BZ[upper.tri(net)], sig, Inf))) +
      sum(log(dst(net[lower.tri(net)], BZ[lower.tri(net)], sig, Inf)))
  }

  log.li <- rbind(theta) %*%
    calculate_h_statistics(
      GERGM_Object,
      GERGM_Object@statistic_auxiliary_data,
      all_weights_are_one = FALSE,
      calculate_all_statistics = FALSE,
      use_constrained_network = TRUE,
      network = net2) +
    last_term
  return(as.numeric(log.li))
}

# maximum pseudo-likelihood estimates
mple <- function(net,
                 statistics,
                 directed,
                 alphas,
                 together,
                 weighted_MPLE,
                 verbose = TRUE) {
  xy <- net2xy(net, statistics, directed, alphas, together)
  x <- xy$x
  y <- xy$y
  # why do we select this initialization
  est <- coef(lm(y ~ x - 1))
  ests <- NULL
  if (verbose) {
      ests <- optim(par = est,
                    pl,
                    y = y,
                    x = x,
                    method = "BFGS",
                    hessian = TRUE,
                    control = list(fnscale = -1, trace = 6))
  } else {
      ests <- optim(par = est,
                    pl,
                    y = y,
                    x = x,
                    method = "BFGS",
                    hessian = TRUE,
                    control = list(fnscale = -1, trace = 0))
  }
  return(ests)
}

# Log likelihood function calculations

# pseudolikelihood given theta#
pl <- function(theta, y, x) {
  return(sum(log(dtexp(y, x %*% theta))))
}

# The conditional density of each weight from a sample
dtexp <- function(x, lambda) {
  den <- numeric(length(x))
  inds <- which(lambda != 0)
  den[inds] <- exp(x[inds] * lambda[inds]) /
    (1 / lambda[inds] * (exp(lambda[inds]) - 1))
  den[which(lambda == 0)] <- 1
  return(den)
}

# The conditional density of each weight from a sample but incorporating alpha
# downweighting. In this version we have to numerically integrate over [0,1] for
# the edge weights.


mple_weighted <- function(GERGM_Object,
                 statistics,
                 possible.stats,
                 verbose = TRUE,
                 prev_ests = NULL) {
  net <- GERGM_Object@network
  num_nodes <- nrow(net)
  triples = t(combn(1:num_nodes, 3))
  pairs <- t(combn(1:num_nodes, 2))

  # we are removing this for now:
  # if (is.null(prev_ests)) {
  #   xy <- net2xy(net,
  #                statistics,
  #                directed = GERGM_Object@directed_network,
  #                alphas = rep(1, length(possible.stats)),
  #                together = 1)
  #   x <- xy$x
  #   y <- xy$y
  #   # why do we select this initialization
  #   est <- coef(lm(y ~ x - 1))
  # } else {
  #   est <- prev_ests
  # }

  est <- prev_ests
  ests <- NULL
  if (verbose) {
    ests <- optim(par = est,
                  extended_fast_pl_weighted,
                  triples = triples,
                  pairs = pairs,
                  GERGM_Object = GERGM_Object,
                  method = "BFGS",
                  hessian = TRUE,
                  control = list(fnscale = -1, trace = 6))
  } else {
    ests <- optim(par = est,
                  extended_fast_pl_weighted,
                  triples = triples,
                  pairs = pairs,
                  GERGM_Object = GERGM_Object,
                  method = "BFGS",
                  hessian = TRUE,
                  control = list(fnscale = -1, trace = 0))
  }
  return(ests)
}

# this version is now deprecated
fast_pl_weighted <- function(theta,
                             triples,
                             pairs,
                             GERGM_Object,
                             lower = 0,
                             upper = 1,
                             steps = 150){

  cat("Weighted MPLE Theta = ",theta,"\n")
  current_network <- GERGM_Object@bounded.network
  number_of_nodes <- nrow(current_network)
  stat.indx <- which(GERGM_Object@stats_to_use > 0)
  #cat("stat.idx",stat.indx,"\n" )
  full_thetas <- rep(0, length(GERGM_Object@stats_to_use))
  for (i in 1:length(theta)) {
    full_thetas[stat.indx[i]] <- theta[i]
  }
  dw <- as.numeric(GERGM_Object@downweight_statistics_together)

  integration_interval <- seq(from = lower,
                              to = upper,
                              length.out = steps)

  objective <- weighted_mple_objective(number_of_nodes,
                            GERGM_Object@stats_to_use,
                            current_network,
                            full_thetas,
                            triples - 1,
                            pairs - 1,
                            GERGM_Object@weights,
                            dw,
                            integration_interval,
                            GERGM_Object@parallel)
  cat("Calculation complete, objective is:",objective,"\n\n")
  return(objective)
}

# current version that works with all of our flexible new statistics
extended_fast_pl_weighted <- function(theta,
                             triples,
                             pairs,
                             GERGM_Object,
                             lower = 0,
                             upper = 1,
                             steps = 150){

  cat("Weighted MPLE Theta = ",theta,"\n")
  current_network <- GERGM_Object@bounded.network
  number_of_nodes <- nrow(current_network)


  dw <- as.numeric(GERGM_Object@downweight_statistics_together)

  integration_interval <- seq(from = lower,
                              to = upper,
                              length.out = steps)

  sad <- GERGM_Object@statistic_auxiliary_data

  # here we are only going to work with our actual thetas/statistics, not the full suite.
  selected_rows_matrix <- sad$specified_selected_rows_matrix
  rows_to_use <- sad$specified_rows_to_use
  base_statistics_to_save <- sad$specified_base_statistics_to_save
  base_statistic_alphas <- sad$specified_base_statistic_alphas
  num_non_base_statistics <- sum(GERGM_Object@non_base_statistic_indicator)

  objective <- extended_weighted_mple_objective(
    number_of_nodes,
    GERGM_Object@stats_to_use - 1,
    current_network,
    triples - 1,
    pairs - 1,
    selected_rows_matrix - 1,
    rows_to_use - 1,
    base_statistics_to_save - 1,
    base_statistic_alphas,
    num_non_base_statistics,
    GERGM_Object@non_base_statistic_indicator,
    theta,
    GERGM_Object@weights,
    dw,
    integration_interval,
    GERGM_Object@parallel)
  cat("Calculation complete, objective is:",objective,"\n\n")
  return(objective)
}


# Convert an observed network to edge weight vectors x and y
net2xy <- function(net, statistics, directed, alphas, together) {
  y <- NULL
  x <- NULL
  nodes <- nrow(net)
  if (directed == TRUE) {
    for (i in 1:nodes) {
      for (j in (1:nodes)[-i]) {
        y <- c(y, net[i, j])
        x <- rbind(x, dh(net, statistics, i, j, alphas, together))
      }
    }
  }
  if (directed == FALSE) {
    for (i in 1:nodes) {
      for (j in (1:nodes)[-i]) {
        y <- c(y, net[i, j])
        x <- rbind(x, dh(net, statistics, i, j, alphas, together))
      }
    }
  }
  return(list(y = y, x = x))
}

# why did we stop doing this in the undirected case?
#     for (i in 2:nodes) {
#       for (j in (1:(i - 1))) {
#         y <- c(y, net[i, j])
#         x <- rbind(x, dh(net, statistics, i, j))
#       }
#     }

# ------------------------------------------------------------
## Functions for correlation matrices
#1: jacobian of transformation of correlation matrices to the [0,1] space
jacobian <- function(partials){
  corrs.1 <- diag(partials[-nrow(partials), -1])
  d <- nrow(partials)
  prod.1 <- prod((1 - corrs.1^2)^(d-2))
  prod.2 <- 1
  for(k in 2 : (d - 2)){
    for(i in 1 : (d - k)){
      prod.2 = prod.2*(1-(partials[i,i+k])^2)^(d-1-k)
    }
  }
  result <- 2*((prod.1^(d-2))*prod.2)^(0.5)
  return(result)
}

#pseudo-likelihood for the correlation matrix
pl.corr <- function(theta, y, x, Jacobian){
  return(sum(log(dtexp(y, x %*% theta))) + log(Jacobian))
}

#MPLE for correlation matrices
mple.corr <- function(net,
                      bounded.net,
                      statistics,
                      alphas,
                      together,
                      directed = FALSE,
                      verbose = TRUE){
  xy.full <- net2xy(net,
                    statistics,
                    directed,
                    alphas,
                    together)
  x <- xy.full$x #x's are the change statistics associated with the unbounded network
  xy.bounded <- net2xy(bounded.net,
                       statistics,
                       directed,
                       alphas,
                       together)
  y <- xy.bounded$y #y's are the edge weights from the bounded [0,1] network
  J <- jacobian(bounded.net)
  est <- coef(lm(y ~ x - 1))
  ests <- NULL
  if(verbose){
    ests <- optim(par = est,
                  pl.corr,
                  y = y,
                  x = x,
                  Jacobian = J,
                  method = "BFGS",
                  hessian = TRUE,
                  control = list(fnscale = -1, trace = 6))
  }else{
    ests <- optim(par = est,
                  pl.corr,
                  y = y,
                  x = x,
                  Jacobian = J,
                  method = "BFGS",
                  hessian = TRUE,
                  control = list(fnscale = -1, trace = 0))
  }
  return(ests)
}

