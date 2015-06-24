# Directional derivative statistics next These are needed for the Gibbs sampler

# din2star
din2star <- function(i, j, net) {
  nodes <- nrow(net)
  others <- (1:nodes)[-c(i, j)]
  sum(net[cbind(others, j)])
}

# dout2star
dout2star <- function(i, j, net) {
  nodes <- nrow(net)
  others <- (1:nodes)[-c(i, j)]
  sum(net[cbind(i, others)])
}

# dedgeweight
dedgeweight = function(i, j) {
  1
}

# drecip
drecip <- function(i, j, net) {
  net[j, i]
}

# dctriads
dctriads <- function(i, j, net) {
  nodes <- nrow(net)
  others <- (1:nodes)[-c(i, j)]
  triples <- cbind(i, j, others)
  sum(net[triples[, c(2, 3)]] * net[triples[, c(3, 1)]])
}

# dttriads
dttriads <- function(i, j, net) {
  nodes <- nrow(net)
  others <- (1:nodes)[-c(i, j)]
  triples <- cbind(i, j, others)
  t2 <- sum(net[triples[, c(2, 3)]] * net[triples[, c(1, 3)]])
  t3 <- sum(net[triples[, c(3, 2)]] * net[triples[, c(3, 1)]])
  t4 <- sum(net[triples[, c(3, 2)]] * net[triples[, c(1, 3)]])
  return(t2 + t3 + t4)
}



# dh function weight w_{i,j} will be conditioned upon
# Calculate the marginal change in the network
dh <- function(net, statistics, i, j) {
  temp <- c(dout2star(i, j, net), din2star(i, j, net), dctriads(i, j, net),
            drecip(i, j, net), dttriads(i, j, net), dedgeweight(i, j))
  value <- temp[statistics > 0]
  return(value)
}





# Functions needed for likelihood function calculations

dgam01 <- function(x, exbz, exlam) {
  x <- x + 0.1
  return(dgamma(x, shape = exbz / exlam, scale = exlam) /
           (1 - pgamma(0.1, shape = exbz / exlam, scale = exlam)))
}

pgam01 <- function(x, exbz, exlam) {
  x <- x + 0.1
  return((pgamma(x, shape = exbz / exlam, scale = exlam) -
            pgamma(0.1, shape = exbz / exlam, scale = exlam)) /
           (1 - pgamma(0.1, shape = exbz / exlam, scale = exlam)))
}

qgam01 <- function(p, exbz, exlam) {
  x <- qgamma(pgamma(0.1, shape = exbz / exlam, scale = exlam) +
                p * (1 - pgamma(0.1, shape = exbz / exlam, scale = exlam)),
              shape = exbz / exlam, scale = exlam)
  return(x - 0.1)
}

dst <- function(x, m, sig, df) {
  return(1 / sig * dt((x - m) / sig, df))
}

pst <- function(x, m, sig, df) {
  return(pt((x - m) / sig, df))
}

qst <- function(u, m, sig, df) {
  return(qt(u, df) * sig + m)
}

indeg <- function(net) {
  id <- numeric(nrow(net))
  for (i in 1:length(id)) {
    id[i] <- sum(net[, i])
  }
  return(id)
}

outdeg <- function(net) {
  od <- numeric(nrow(net))
  for (i in 1:length(od)) {
    od[i] <- sum(net[i, ])
  }
  return(od)
}

# Convert an observed network to edge weight vectors x and y

net2xy <- function(net, statistics, directed) {
  y <- NULL
  x <- NULL
  nodes <- nrow(net)
  if (directed == TRUE) {
    for (i in 1:nodes) {
      for (j in (1:nodes)[-i]) {
        y <- c(y, net[i, j])
        x <- rbind(x, dh(net, statistics, i, j))
      }
    }

  }
  if (directed == FALSE) {
    for (i in 2:nodes) {
      for (j in (1:(i - 1))) {
        y <- c(y, net[i, j])
        x <- rbind(x, dh(net, statistics, i, j))
      }
    }

  }
  return(list(y = y, x = x))
}

# Draw a random value either uniform or according to density
rtexp <- function(n, lambda) {
  # lambda is the scalar-valued parameter n is the number of draws
  u <- runif(n)
  if (lambda != 0) {
    temp = exp(lambda)
    temp[is.infinite(temp)] = 1e+32
    x <- log(1 + u * (temp - 1)) / lambda
    x
  } else u
}

# The conditional density of each weight from a sample
dtexp <- function(x, lambda) {
  den <- numeric(length(x))
  den[which(lambda != 0)] <- exp(x[which(lambda != 0)] *
                                   lambda[which(lambda != 0)]) /
    (1 / lambda[which(lambda != 0)] * (exp(lambda[which(lambda != 0)]) - 1))
  den[which(lambda == 0)] <- 1
  return(den)
}

# Function to generate dispersed unit interval
rdisp <- function(n) {
  pnorm(abs(rnorm(n) * 6))
}

# Log likelihood function calculations

# pseudolikelihood given theta#
pl <- function(theta, y, x) {
  return(sum(log(dtexp(y, x %*% theta))))
}

# add to console output field in GERGM_Object
#GERGM_Object <- store_console_output(GERGM_Object,addition)
store_console_output <- function(GERGM_Object,addition){
  if("list" %in% class(addition)){
    addition <- as.character(unlist(addition))
  }
  if(is.null(GERGM_Object@console_output)){
    GERGM_Object@console_output <- addition
  }else{
    GERGM_Object@console_output <- c(GERGM_Object@console_output,addition)
  }
  return(GERGM_Object)
}

