# Statistics out2stars
out2star <- function(net, triples, alpha = 1, together = together) {
  if(together == 0){
    st1 <- sum(net[triples[, c(1, 2)]]^(alpha) * net[triples[, c(1, 3)]]^(alpha))
    st2 <- sum(net[triples[, c(2, 1)]]^(alpha) * net[triples[, c(2, 3)]]^(alpha))
    st3 <- sum(net[triples[, c(3, 1)]]^(alpha) * net[triples[, c(3, 2)]]^(alpha))
    return(st1 + st2 + st3)
  }
  if(together == 1){
    st1 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(1, 3)]])
    st2 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(2, 3)]])
    st3 <- sum(net[triples[, c(3, 1)]] * net[triples[, c(3, 2)]])
    return((st1 + st2 + st3)^alpha)
  }
}

# in2stars
in2star <- function(net, triples, alpha = 1, together = together) {
  if(together == 0){
    st1 <- sum(net[triples[, c(3, 1)]]^(alpha) * net[triples[, c(2, 1)]]^(alpha))
    st2 <- sum(net[triples[, c(3, 2)]]^(alpha) * net[triples[, c(1, 2)]]^(alpha))
    st3 <- sum(net[triples[, c(1, 3)]]^(alpha) * net[triples[, c(2, 3)]]^(alpha))
    return(st1 + st2 + st3)
  }
  if(together == 1){
    st1 <- sum(net[triples[, c(3, 1)]] * net[triples[, c(2, 1)]])
    st2 <- sum(net[triples[, c(3, 2)]] * net[triples[, c(1, 2)]])
    st3 <- sum(net[triples[, c(1, 3)]] * net[triples[, c(2, 3)]])
    return((st1 + st2 + st3)^alpha)
  }
}

# transitive triads
ttriads <- function(net, triples, alpha = 1, together) {
  if(together == 0){
    t2 <- sum(net[triples[, c(1, 2)]]^(alpha) * net[triples[, c(2, 3)]]^(alpha) *
                net[triples[, c(1, 3)]]^(alpha))
    t3 <- sum(net[triples[, c(1, 2)]]^(alpha) * net[triples[, c(3, 2)]]^(alpha) *
                net[triples[, c(3, 1)]]^(alpha))
    t4 <- sum(net[triples[, c(1, 2)]]^(alpha) * net[triples[, c(3, 2)]]^(alpha) *
                net[triples[, c(1, 3)]]^(alpha))
    t5 <- sum(net[triples[, c(2, 1)]]^(alpha) * net[triples[, c(2, 3)]]^(alpha) *
                net[triples[, c(3, 1)]]^(alpha))
    t6 <- sum(net[triples[, c(2, 1)]]^(alpha) * net[triples[, c(2, 3)]]^(alpha) *
                net[triples[, c(1, 3)]]^(alpha))
    t7 <- sum(net[triples[, c(2, 1)]]^(alpha) * net[triples[, c(3, 2)]]^(alpha) *
                net[triples[, c(3, 1)]]^(alpha))
    return(t2 + t3 + t4 + t5 + t6 + t7)
  }
  if(together == 1){
    t2 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(2, 3)]] *
                net[triples[, c(1, 3)]])
    t3 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(3, 2)]] *
                net[triples[, c(3, 1)]])
    t4 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(3, 2)]] *
                net[triples[, c(1, 3)]])
    t5 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(2, 3)]] *
                net[triples[, c(3, 1)]])
    t6 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(2, 3)]] *
                net[triples[, c(1, 3)]])
    t7 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(3, 2)]] *
                net[triples[, c(3, 1)]])
    return((t2 + t3 + t4 + t5 + t6 + t7)^alpha)
  }
}

# cyclic triads
ctriads <- function(net, triples, alpha = 1, together) {
  if(together == 0){
    t1 <- sum(net[triples[, c(1, 2)]]^(alpha) * net[triples[, c(2, 3)]]^(alpha) *
                net[triples[, c(3, 1)]]^(alpha))
    t8 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(3, 2)]] *
                net[triples[, c(1, 3)]]^(alpha))
    return(t1 + t8)
  }
  if(together == 1){
    t1 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(2, 3)]] *
                net[triples[, c(3, 1)]])
    t8 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(3, 2)]] *
                net[triples[, c(1, 3)]])
    return((t1 + t8)^alpha)
  }
}

# reciprocity
recip <- function(net, alpha = 1, together) {
  pairs <- t(combn(1:nrow(net), 2))
  if(together == 0){
    return(sum(net[pairs]^(alpha) * net[pairs[, c(2, 1)]]^(alpha)))
  }
  if(together == 1){
    return(sum(net[pairs] * net[pairs[, c(2, 1)]])^(alpha))
  }
}

# edgeweight
edgeweight <- function(net, alpha = 1, together) {
  pairs <- t(combn(1:nrow(net), 2))
  if(together == 0){
    return(sum((net[pairs]^(alpha) + net[pairs[, c(2, 1)]])^(alpha)))
  }
  if(together == 1){
    return(sum((net[pairs] + net[pairs[, c(2, 1)]]))^(alpha))
  }
}

# Calculate the statistics of a formula object
h <- function(formula, possible.stats, alpha = NULL, theta = NULL, together = together) {
  res1 <- Parse_Formula_Object(formula = formula, possible.stats, alpha = alpha, theta = theta)
  net <- res1$net
  alphas <- res1$alphas
  num.nodes <- nrow(net)
  triples = t(combn(1:num.nodes, 3))
  statistics <- res1$statistics
  temp <- c(out2star(net, triples, alphas[1], together),
            in2star(net, triples, alphas[2], together),
            ctriads(net, triples, alphas[3], together),
            recip(net, alphas[4], together),
            ttriads(net, triples, alphas[5], together),
            edgeweight(net, alphas[6], together))
  value <- temp[statistics > 0]
  result <- rbind(round(value, 3), round(alphas[statistics > 0], 3))
  colnames(result) <- possible.stats[statistics > 0]
  rownames(result) <- c("value", "alpha")
  return(result)
}

# A second version of the h function used for calculation in estimation
# This function calculates the network statistics associated with net
h2 <- function(net, triples, statistics, alphas = rep(1,6), together = 1) {
  temp = c(out2star(net, triples, alphas[1], together),
           in2star(net, triples, alphas[2], together),
           ctriads(net, triples, alphas[3], together),
           recip(net, alphas[4], together),
           ttriads(net, triples, alphas[5], together),
           edgeweight(net, alphas[6], together))
  return(temp[which(statistics > 0)])
  #return(temp)
}
