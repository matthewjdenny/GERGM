# maximum pseudo-likelihood estimates
ex_mple <- function(GERGM_Object,
                    verbose = TRUE) {
  xy <- ex_net2xy(GERGM_Object)
  x <- xy$x
  y <- xy$y
  # why do we select this initialization
  est <- coef(lm(y ~ x - 1))
  ests <- NULL
  if (verbose) {
    ests <- optim(par = est,
                  ex_pl,
                  y = y,
                  x = x,
                  method = "BFGS",
                  hessian = TRUE,
                  control = list(fnscale = -1, trace = 6))
  } else {
    ests <- optim(par = est,
                  ex_pl,
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
ex_pl <- function(theta, y, x) {
  return(sum(log(ex_dtexp(y, x %*% theta))))
}

# The conditional density of each weight from a sample
ex_dtexp <- function(x, lambda) {
  den <- numeric(length(x))
  inds <- which(lambda != 0)
  den[inds] <- exp(x[inds] * lambda[inds]) /
    (1 / lambda[inds] * (exp(lambda[inds]) - 1))
  den[which(lambda == 0)] <- 1
  return(den)
}


# Convert an observed network to edge weight vectors x and y
ex_net2xy <- function(GERGM_Object) {
  net <- GERGM_Object@bounded.network
  y <- NULL
  x <- NULL
  nodes <- GERGM_Object@num_nodes
  for (i in 1:nodes) {
    for (j in (1:nodes)[-i]) {
      y <- c(y, net[i, j])
      x <- rbind(x, ex_dh(GERGM_Object, i, j))
    }
  }
  return(list(y = y, x = x))
}


# din2star
ex_din2star <- function(i, j, net, node_set) {
  val <- 0
  if (i %in% node_set & j %in% node_set) {
    others <- node_set[-c(i, j)]
    val <- sum(net[cbind(others, j)])
  }
  return(val)
}

# dout2star
ex_dout2star <- function(i, j, net, node_set) {
  val <- 0
  if (i %in% node_set & j %in% node_set) {
    others <- node_set[-c(i, j)]
    val <- sum((net[cbind(i, others)]))
  }
  return(val)
}

# dedgeweight
ex_dedgeweight <- function(i, j, net, node_set) {
  val <- 0
  if (i %in% node_set & j %in% node_set) {
    val <- 1
  }
  return(val)
}

# drecip
ex_drecip <- function(i, j, net, node_set) {
  val <- 0
  if (i %in% node_set & j %in% node_set) {
    val <- net[j, i]
  }
  return(val)
}

# dctriads
ex_dctriads <- function(i, j, net, node_set) {
  val <- 0
  if (i %in% node_set & j %in% node_set) {
    others <- node_set[-c(i, j)]
    triples <- cbind(i, j, others)
    val <- sum(net[triples[, c(2, 3)]] * net[triples[, c(3, 1)]])
  }
  return(val)
}

# dttriads
ex_dttriads <- function(i, j, net, node_set) {
  val <- 0
  if (i %in% node_set & j %in% node_set) {
    others <- node_set[-c(i, j)]
    triples <- cbind(i, j, others)
    t2 <- sum(net[triples[, c(2, 3)]] * net[triples[, c(1, 3)]])
    t3 <- sum(net[triples[, c(3, 2)]] * net[triples[, c(3, 1)]])
    t4 <- sum(net[triples[, c(3, 2)]] * net[triples[, c(1, 3)]])
    val <- t2 + t3 + t4
  }
  return(val)
}

# dh function weight w_{i,j} will be conditioned upon
# Calculate the marginal change in the network
ex_dh <- function(GERGM_Object, i, j, network = NULL) {
  # find out what stats we are using
  stats_to_use <- GERGM_Object@stats_to_use
  value <- rep(0, length(stats_to_use))

  # allow the user to pass in a network
  if (!is.null(network)) {
    net <- network
  } else {
    net <- GERGM_Object@bounded.network
  }

  # now loop through them to calculate d statistics
  for (k in 1:length(stats_to_use)) {
    # determine the node set
    node_set <- GERGM_Object@endogenous_statistic_node_sets[[k]]
    # if there was no reduced nodeset provided, then use all nodes
    if (length(node_set) == 1) {
      node_set <- 1:GERGM_Object@num_nodes
    }
    # now get the value depending on which statistic we are dealing with
    if (stats_to_use[k] == 1) {
      value[k] <- ex_dout2star(i, j, net, node_set)
    } else if (stats_to_use[k] == 2) {
      value[k] <- ex_din2star(i, j, net, node_set)
    } else if (stats_to_use[k] == 3) {
      value[k] <- ex_dctriads(i, j, net, node_set)
    } else if (stats_to_use[k] == 4) {
      value[k] <- ex_drecip(i, j, net, node_set)
    } else if (stats_to_use[k] == 5) {
      value[k] <- ex_dttriads(i, j, net, node_set)
    } else if (stats_to_use[k] == 6) {
      value[k] <- ex_dedgeweight(i, j, net, node_set)
    } else {
      stop("MPLE is not currently implemented for this statistic.")
    }
  }
  return(value)
}
