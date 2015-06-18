# Last update: April 12th, 2015

# an S4 class for gergm objects
setClass(Class = "gergm", 
         representation = representation(
           network = "matrix",
           bounded.network = "matrix",
           formula = "formula",
           stats = "matrix",
           theta.coef = "data.frame", 
           lambda.coef = "data.frame",
           weights = "numeric"
         ), 
         validity = function(object) {
           if (!"matrix" %in% class(object@network) & is.null(object@network) 
               == FALSE) {
             stop("'network' must be a 'matrix' object or 'NULL'.")
           }
           if (!"matrix" %in% class(object@bounded.network)) {
             stop("'bounded.network' must be a 'matrix' object.")
           }
           if (!is.data.frame(object@theta.coef)) {
             stop("'theta.coef' must be a 'data.frame' object.")
           }
           if (!"formula" %in% class(object@formula)) {
             stop("'formula' is not a 'formula' object.")
           }
           if (!is.data.frame(object@lambda.coef) & is.null(object@lambda.coef) 
               == FALSE) {
             stop("'lambda.coef' must be a 'data.frame' object or 'NULL'.")
           }
           if (!"matrix" %in% class(object@stats)) {
             stop("'stats' must be a 'matrix' object.")
           }
           if (!is.numeric(object@weights) & is.null(object@weights) 
               == FALSE) {
             stop("'weights' must be a 'data.frame' object or 'NULL'.")
           }
           return(TRUE)
         }
)

# define coef for pretty output of gergm object
setMethod(f = "coef", signature = "gergm", definition = function(object, ...) {
  return(list(Theta  = object@theta.coef, Lambda = object@lambda.coef))
}
)

# define 'show' to get a pretty output of the gergm
setMethod(f = "show", signature = "gergm", definition = function(object){
  message("Theta:")
  print(object@theta.coef)
  message("Lambda:")
  print(object@lambda.coef)
  message("Weights:")
  print(object@weights)
  message("Network Statistics:")
  print(object@stats)
}
)

# gergm object constructor
gergm.object <- function(network, bounded.network = network, formula, 
                         thetas, lambda, alpha, together = together) {
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


# create a gergm from a formula object
getgergm <- function(object, theta.coef, together = 1, weights = NULL, 
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

# Simulate a gergm
simulate.gergm <- function(object, nsim, seed = NULL,
                           method = c("Gibbs", "Metropolis"), MCMC.burnin = 100,  
                           weights = object@weights,
                           coef = object@theta.coef, formula = object@formula,
                           thin = 1, shape.parameter = 1, together = 0,
                           seed1 = 10000) {
  # object: an object of class "gergm"
  method <- method[1]
  
  theta <- as.numeric(coef[1,])  # obtain the estimated thetas
  
  sample_every <- floor(1/thin)
  
  net <- object@bounded.network  # obtain the bounded network object
  alpha <- weights
  rhs <- names(object@theta.coef)[abs(theta) > 0]
  rhs <- paste(rhs, collapse = "+")
  formula.new <- as.formula(paste0("net ~ ", rhs))
  
  res1 <- weighted.ergm.data(formula.new, theta = theta[abs(theta) > 0], 
                             alpha = alpha[abs(theta) > 0])
  statistics <- res1$statistics
  
  alphas <- res1$alphas
  thetas <- res1$thetas
  num.nodes <- nrow(net)
  triples <- t(combn(1:num.nodes, 3))
  pairs <- t(combn(1:num.nodes, 2))
  
  # Initial Network is a graph with uniform random edgeweights
  #initial_network = matrix(runif(num.nodes^2, 1e-05, 0.99999), nrow = num.nodes, 
  #ncol = num.nodes)
  initial_network <- net
  
  # Gibbs Simulation
  if (method == "Gibbs") {
    nets <- gergm.gibbs.sampler(formula.new, thetas[statistics > 0], 
                                MCMC.burnin = MCMC.burnin, num.draws = nsim, 
                                thin = thin, start = NULL, num.nodes = num.nodes,
                                directed = TRUE)
    # Calculate the network statistics over all of the simulated networks
    h.statistics <- t(apply(nets, 3, h2, triples = triples,
                            statistics = rep(1, 6), alphas = rep(1, 6),
                            together = together))
    acceptance.rate <- NULL
  }
  
  # Metropolis Hastings Simulation
  if (method == "Metropolis") {
    #fix alphas for MH sims
    alphas[which(is.na(alphas))] = 1
    alphas[which(alphas == 0)] = 1
    #print to check what parameters are going in to the MH sampler
    #print(Sys.time())
    #print("Statistics to use")
    #print(statistics)
    print("Thetas")
    print(thetas)
    #print("Alphas")
    #print(alphas)
    #print("Take sample every:")
    #print(sample_every)
    #print("Shape parameter:")
    #print(shape.parameter)
    #print("Burnin")
    #print(MCMC.burnin)
    #print("Number of Nodes")
    #print(num.nodes)
    #print("together")
    #print(together)
    #print("seed")
    #print(seed1)
    
    samples <- MH(
      number_of_iterations = nsim + MCMC.burnin, 
      shape_parameter = shape.parameter, number_of_nodes = num.nodes, 
      statistics_to_use = statistics, initial_network = initial_network, 
      take_sample_every = sample_every, thetas = thetas, 
      triples = triples - 1, pairs = pairs - 1, alphas = alphas, 
      together = together, seed = seed1)
    # keep only the networks after the burnin
    start <- floor(MCMC.burnin/sample_every) + 1
    end <- length(samples[[3]][,1])
    nets <- samples[[2]][, , start:end]
    # Calculate the network statistics over all of the simulated networks
    # Note: these statistics will be the adjusted statistics (for use in the
    # MCMCMLE procedure)
    
    h.statistics <- samples[[3]][start:end,]
    
    
    acceptance.rate <- mean(samples[[1]])
    cat("acceptance.rate", acceptance.rate, "\n")
    
    # hack to make sure we do not have severe autocorrelation in our sample
    #densities <- h.statistics[,1]
    #cat("number of raw samples:",length(densities),"\n")
    #ar1 <- cor(densities[2:length(densities)],densities[1:(length(densities)-1)])
    #print(ar1)
    
    #         thin <- 1
    #         while(abs(ar1) > .3){
    #             thin = thin +1
    #             thinSeq <- round(seq(1,length(densities),by=thin))
    #             thinDens <- densities[thinSeq]
    #             ar1 <- cor(thinDens[2:length(thinDens)],thinDens[1:(length(thinDens)-1)])
    #             #print(paste(ar1,thin))
    #         }
    #         print(thin)
    #         print(ar1)
    #thin = 200
    
    #thinSeq <- round(seq(1,length(h.statistics[,1]),by=thin))
    
    #h.statistics <- h.statistics[thinSeq,]
    
    
    #h.statistics <- t(apply(nets, 3, h2, triples = triples,
    #                        statistics = rep(1, 6), together = together))
    #hack to add measurement error
    #h.statistics <- h.statistics + runif(6*nsim, min = -100, max = 100)
    
  }
  h.statistics = data.frame(out2stars = h.statistics[, 1], 
                            in2stars = h.statistics[, 2], 
                            ctriads = h.statistics[, 3], 
                            recip = h.statistics[, 4], 
                            ttriads = h.statistics[, 5], 
                            edgeweight = h.statistics[, 6])
  #possible.stats <- c("out2stars", "in2stars", "ctriads", "recip",
  #  "ttriads", "edgeweight")
  #colnames(h.statistics) <- possible.stats
  #h.statistics <- as.data.frame(h.statistics)
  return(list(Networks = nets, Statistics = h.statistics,
              Acceptance.rate = acceptance.rate))
}


# Statistics out2stars
out2star <- function(net, triples, alpha = 1, together = together) {
  if(together == 0){
    st1 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(1, 3)]]^(alpha))
    st2 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(2, 3)]]^(alpha))
    st3 <- sum(net[triples[, c(3, 1)]] * net[triples[, c(3, 2)]]^(alpha))
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
    st1 <- sum(net[triples[, c(3, 1)]] * net[triples[, c(2, 1)]]^(alpha))
    st2 <- sum(net[triples[, c(3, 2)]] * net[triples[, c(1, 2)]]^(alpha))
    st3 <- sum(net[triples[, c(1, 3)]] * net[triples[, c(2, 3)]]^(alpha))
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
    t2 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(2, 3)]] * 
                net[triples[, c(1, 3)]]^(alpha))
    t3 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(3, 2)]] * 
                net[triples[, c(3, 1)]]^(alpha))
    t4 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(3, 2)]] * 
                net[triples[, c(1, 3)]]^(alpha))
    t5 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(2, 3)]] * 
                net[triples[, c(3, 1)]]^(alpha))
    t6 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(2, 3)]] * 
                net[triples[, c(1, 3)]]^(alpha))
    t7 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(3, 2)]] * 
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
    t1 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(2, 3)]] * 
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
    return(sum(net[pairs] * net[pairs[, c(2, 1)]]^(alpha)))
  }
  if(together == 1){
    return(sum(net[pairs] * net[pairs[, c(2, 1)]])^(alpha))
  } 
}

# edgeweight
edgeweight = function(net, alpha = 1, together) {
  pairs <- t(combn(1:nrow(net), 2))
  if(together == 0){
    return(sum((net[pairs]^(alpha) + net[pairs[, c(2, 1)]])^(alpha)))
  }
  if(together == 1){
    return(sum((net[pairs] + net[pairs[, c(2, 1)]]))^(alpha))
  }
}


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

# Calculate the statistics of a formula object
h <- function(formula, alpha = NULL, theta = NULL, together = together) {
  res1 <- weighted.ergm.data(formula = formula, alpha = alpha, theta = theta)
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
  possible.stats <- c("out2star", "in2star", "ctriads", "recip", 
                      "ttriads", "edgeweight")
  result <- rbind(round(value, 3), round(alphas[statistics > 0], 3))
  colnames(result) <- possible.stats[statistics > 0]
  rownames(result) <- c("value", "alpha")
  return(result)
}

# dh function weight w_{i,j} will be conditioned upon
# Calculate the marginal change in the network
dh <- function(net, statistics, i, j) {
  temp <- c(dout2star(i, j, net), din2star(i, j, net), dctriads(i, j, net), 
            drecip(i, j, net), dttriads(i, j, net), dedgeweight(i, j))
  value <- temp[statistics > 0]
  return(value)
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

# convert a formula object to statistics and network of interest
weighted.ergm.data <- function(formula, theta = NULL, alpha = NULL) {
  # parse the formula
  if (class(formula) != "formula") {
    stop("'formula' must be a formula object.")
  }
  lhs <- deparse(formula[[2]])  # name of the response variable
  net <- eval(parse(text = lhs))  # get the actual response data
  
  # check that the response is a matrix
  if (class(net) != "matrix") {
    stop("Response must be a matrix object.")
  }
  
  rhs <- paste0(deparse(formula[[3]]), collapse = "")  # rhs of formula
  rhs <- gsub("\\s+", " ", rhs)  # get rid of redundant spaces
  rhs <- strsplit(rhs, " \\+ ")[[1]]  # parse separate formula elements
  
  possible.stats = c("out2star", "in2star", "ctriads", "recip",
                     "ttriads", "edgeweight")
  
  # check that the names of the statistics match those that are possible
  possible <- 1:length(rhs)
  actual <- which(rhs %in% possible.stats)
  if (length(possible) != length(actual)) {
    stop(paste("the specified variable", rhs[setdiff(possible, actual)],
               "is not an available statistic.", sep = " "))
  }
  
  # if alpha is NULL, assume all ones
  if (is.null(alpha) == TRUE) {
    alpha <- rep(1, length(rhs))
  }
  
  # if theta is NULL, assume all ones
  if (is.null(theta) == TRUE) {
    theta <- rep(1, length(rhs))
  }
  # check that alpha and the number of statistics are equal
  if (length(rhs) != length(alpha)) {
    stop("'alpha' must be the same length as the number of statistics")
  }
  
  # check that theta and the number of statistics are equal
  if (length(rhs) != length(theta)) {
    stop("'theta' must be the same length as the number of statistics")
  }
  stat.indx <- which(possible.stats %in% rhs)
  statistics <- rep(0, 6)
  statistics[stat.indx] <- 1
  alphas <- rep(1, 6)
  thetas <- rep(0, 6)
  for (i in 1:length(rhs)) {
    alphas[which(rhs[i] == possible.stats)] <- alpha[i]
    thetas[which(rhs[i] == possible.stats)] <- theta[i]
  }
  return(list(net = net, statistics = statistics, 
              alphas = alphas, thetas = thetas))
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

# Simulation Functions
gergm.gibbs.sampler <- function(formula, theta, MCMC.burnin, num.draws, 
                                thin = 1, start = NULL, num.nodes = NULL, 
                                directed = c(TRUE, FALSE)) {
  # formula specifies which variables to use MCMC.burnin is the number of
  # discarded draws n is the total number of draws to be reported dh is a 
  # function that takes as arguments, i, j, theta, net and returns the partial 
  # of the hamiltonian theta is the vector-valued parameter thin reduces 
  # autocorrelation in the simulations, every thinth iteration is returned start
  # is the initial network, if not supplied, a random uniform nodesXnodes 
  # network is used num.nodes is the number of nodes in the network dir is a 
  # logical indicator of whether the network is directed
  
  # Extract which statistics we use
  res1 <- weighted.ergm.data(formula, theta = theta, alpha = NULL)
  net <- res1$net
  if (is.null(num.nodes) == TRUE) {
    num.nodes <- nrow(net)
  }
  thetas <- res1$thetas
  statistics <- res1$statistics
  sims <- num.draws * thin
  netarray <- array(NA, dim = c(num.nodes, num.nodes, num.draws + 1))
  if (is.null(start)) 
    start <- matrix(rdisp(num.nodes * num.nodes), num.nodes, num.nodes)
  net <- start
  diag(net) = 0
  netarray[, , 1] <- net
  for (t in 1:(num.draws + MCMC.burnin)) {
    if (directed == TRUE) {
      for (i in 1:num.nodes) {
        for (j in (1:num.nodes)[-i]) {
          net[i, j] <- rtexp(1, t(theta) %*% dh(net, statistics, i, j))
        }
      }
    }
    if (directed == FALSE) {
      for (i in 2:num.nodes) {
        for (j in (1:(i - 1))) {
          net[i, j] <- rtexp(1, t(theta) %*% dh(net, statistics, i, j))
        }
      }
    }
    if (t > MCMC.burnin) {
      netarray[, , t - MCMC.burnin] <- net
    }
  }
  return(netarray[, , round(seq(1, sims, length = num.draws))])
}

# Estimation Code

# log likelihood
log.l <- function(thetas, formula, alpha, hsnet, ltheta, together = together) {
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
  
  temp <- h(formula, alpha, theta = theta, together = together)[1, ]
  return(rbind(theta) %*% temp - max(z) - log(sum(exp(z - max(z)))))
  
  #theta <- par[1:ncol(hsnet)]
  #z <- hsnet%*%(theta-ltheta)
  #rbind(theta)%*%h(net,triples)-max(z)-log(sum(exp(z-max(z))))
  #return(rbind(theta) %*% temp - max(z) - log(mean(exp(z - max(z)))))
}

llg <- function(par, formula, alpha, theta, z, together = together) {
  # log likelihood for unbounded network with g function
  res1 <- weighted.ergm.data(formula, theta = theta, alpha = alpha)
  statistics <- res1$statistics
  alphas <- res1$alphas
  #cat("alphas in llg")
  #print(alphas)
  net <- res1$net
  #cat("net in llg")  #Network is correct
  #print(head(net))
  #theta <- res1$thetas
  #possible.stats <- c("out2star", "in2star", "ctriads", "recip",
  # "ttriads", "edgeweight")
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

# The main MCMC maximum likelihood estimation step.  
# Swaps between simulation and estimation.
# This only concerns graphs with weights between [0,1]

mcmcmle <- function(formula.obj, num.draws, mc.num.iterations,
                    tolerance, thin = 1, MCMC.burnin, theta = NULL, 
                    alpha = NULL, directed = c(TRUE, FALSE), 
                    method = c("Gibbs", "Metropolis"), shape.parameter = 1, 
                    take.sample.every = 1, together = 0, seed2 = 10000, gain.factor) {
  directed <- directed[1]
  method <- method[1]
  res1 <- weighted.ergm.data(formula.obj, theta = theta, alpha = alpha)
  statistics <- res1$statistics
  alphas <- res1$alphas
  #cat("alphas in mcmcmle")
  #print(alphas)
  net2 <- res1$net 
  
  if(CUSTOM){
    theta.init <- list()
    theta.init$par <- MYTHETAS
  } else{
    theta.init <- mple(net2, statistics = statistics, directed = directed)
  }
  # print the Theta parameters
  #theta$par <- My_Parameters
  #cat("Initial Thetas: ", theta$par, "\n")
  cat("MPLE Thetas: ", theta.init$par, "\n")
  num.nodes <- nrow(net2)
  triples <- t(combn(1:num.nodes, 3))
  pairs <- t(combn(1:nrow(net2), 2))
  # initialize the network with the observed network
  initial_network <- net2
  # calculate the statistics of the original network
  init.statistics <- h2(net2, triples = triples, statistics = rep(1, 6),
                        alphas = alphas, together = together)
  obs.stats <- h2(net2, triples = triples, statistics = statistics, alphas = alphas, together = together)
  
  cat("Observed statistics", "\n", obs.stats, "\n")
  #################################################################################################
  ##JW: Added 3/29/15. This scales the initial estimates for the MPLE theta specification
  ## This is according to the initialization the Fisher Scoring method for optimization
  alps <- alphas[which(statistics == 1)]
  
  object <- getgergm(formula.obj, theta.coef = theta.init$par, 
                     weights = alps, together = together)
  temp <- simulate.gergm(object, nsim = ceiling(20/thin), method = method, 
                         shape.parameter = shape.parameter, 
                         together = together, thin = thin,
                         MCMC.burnin = MCMC.burnin, seed1 = seed2)
  
  hsn <- temp$Statistics[,which(statistics == 1)]
  
  #Set gain factor. Here, we will set the default of 0.10, but this may not be enough!
  #gain.factor <- 0.01
  #Calculate covariance estimate (to scale initial guess theta.init)
  z.bar <- colSums(hsn) / 20
  cat("z.bar", "\n", z.bar, "\n")
  Cov.est <- 0
  for(i in 1:dim(hsn)[1]){
    Cov.est <- matrix(as.numeric(hsn[i,]), ncol = 1) %*% t(matrix(as.numeric(hsn[i,]), ncol = 1)) + Cov.est
  }
  Cov.est <- (Cov.est / 20) - z.bar%*%t(z.bar)
  #cat("Cov.est", "\n", Cov.est)
  D.inv <- solve(Cov.est)
  #calculate
  theta <- list()
  theta$par <- theta.init$par - gain.factor * D.inv %*% (z.bar - obs.stats) 
  cat("Adjusted initial theta", "\n", theta$par, "\n")
  #################################################################################################  
  
  ## Simulate new networks
  for (i in 1:mc.num.iterations) {
    alps <- alphas[which(statistics == 1)]
    object <- getgergm(formula.obj, theta.coef = theta$par, 
                       weights = alps, together = together)
    temp <- simulate.gergm(object, nsim = num.draws, method = method, 
                           shape.parameter = shape.parameter, 
                           together = together, thin = thin,
                           MCMC.burnin = MCMC.burnin, seed1 = seed2)
    
    
    #just use what gets returned
    #snets <- temp$Networks 
    #hsn <- t(apply(snets, 3, h2, triples = triples, statistics = statistics,
    #               alphas = alphas, together = together))
    hsn <- temp$Statistics[,which(statistics == 1)]
    
    
    #for(i in 1:length(alps)){
    #    hsn[,i] <- hsn[,i] + runif(length(hsn[,i]),min = -(My_Bounds),max = My_Bounds)
    #}
    
    #hsn.tot <- t(apply(snets, 3, h2, triples = triples, statistics = rep(1, 6),
    #                   alphas = alphas, together = together))
    hsn.tot <- temp$Statistics
    #cat("Simulations Done", "\n")
    #calculate t.test p-values for calculating the difference in the means of 
    # the newly simulated data with the original network
    t.out <- t.test(hsn.tot[, 1], mu = init.statistics[1])$p.value
    t.in <- t.test(hsn.tot[, 2], mu = init.statistics[2])$p.value
    t.ctriad <- t.test(hsn.tot[, 3], mu = init.statistics[3])$p.value
    t.recip <- t.test(hsn.tot[, 4], mu = init.statistics[4])$p.value
    t.ttriads <- t.test(hsn.tot[, 5], mu = init.statistics[5])$p.value
    t.edge <- t.test(hsn.tot[, 6], mu = init.statistics[6])$p.value
    
    stats.data <- data.frame(Observed = init.statistics, 
                             Simulated = colMeans(hsn.tot))
    rownames(stats.data) <- c("out2star", "in2star", "ctriads", "recip",
                              "ttriads", "edgeweight")
    print(stats.data)
    #t.stats <- data.frame(out2star = t.out, in2star = t.in, ctriads = t.ctriad,
    #recip = t.recip, ttriads = t.ttriads, 
    #edgeweight = t.edge)
    #cat("t test p-values:", "\n")
    #print(t.stats)
    #cat("alpha before theta optim", "\n")
    #print(alps)
    #cat("theta before optim")
    #print(theta$par)
    theta.new <- optim(par = theta$par, log.l, alpha = alps, 
                       formula = formula.obj, hsnet = hsn, ltheta = as.numeric(theta$par), 
                       together = together, method = "BFGS", hessian = T, 
                       control = list(fnscale = -1, trace = 5))
    print(paste0("Theta = ", theta.new$par))
    theta.std.errors <- 1 / sqrt(abs(diag(theta.new$hessian)))
    # Calculate the p-value based on a z-test of differences
    # The tolerance is the alpha at which differences are significant
    # JW: Changed on April 1st
    #Comment out the old calculation
    #t.stat <- qt(1 - (tolerance / 2), 1000000)
    #bounds <- t.stat * theta.std.errors
    #cat("Bounds \n", bounds, "\n")
    p.value <- rep(0,length(theta$par))
    count <- rep(0, length(theta$par))
    for(i in 1:length(theta$par)){
      #two sided z test
      p.value[i] <- 2*pnorm(-abs((theta.new$par[i] - theta$par[i])/theta.std.errors[i]))
      #abs(theta.new$par[i] - theta$par[i]) > bounds[i]
      #if we reject any of the tests then convergence has not been reached!
      if(p.value[i] < tolerance){count[i] = 1}
    }
    cat("p.values", "\n")
    print(p.value)
    #cat("Differences", "\n")
    #cat(abs(theta.new$par - theta$par))
    #if (max(abs(theta.new$par - theta$par)) < tolerance) {}
    if (sum(count) == 0){
      message("Parameter estimates have converged")
      return(theta.new)
    }
    else{
      cat("\n", "Theta Estimates", theta.new$par, "\n")
    }
    theta = theta.new
  }
  return(theta)
}

# Function for plotting goodness of fit statistics in a boxplot
gof.plot <- function(simulation.obj, gergm.obj){
  # simulation.obj is a list object resulting from simulations
  # gergm.obj is an object of class "gergm"
  boxplot(simulation.obj$Statistics, medcol = "red")
  boxplot(rbind(gergm.obj@stats[2, ], gergm.obj@stats[2, ]), add =T, 
          medcol="blue", names = F)
}

gof.compare.plot <- function(MH.obj, Gibbs.obj, gergm.obj){
  par(mfrow = c(2,3))
  violins(data.frame(Gibbs.obj$Statistics[,1], MH.obj$Statistics[,1]), names = c("Gibbs", "M-H"), col = c("blue", "green"), connectcol = "transparent", main = "Out-2-Stars")
  abline(h = gergm.obj@stats[2,1], col = "red", lty = 2, lwd = 2)
  
  violins(data.frame(Gibbs.obj$Statistics[,2], MH.obj$Statistics[,2]), names = c("Gibbs", "M-H"), col = c("blue", "green"), connectcol = "transparent", main = "In-2-Stars")
  abline(h = gergm.obj@stats[2,2], col = "red", lty = 2, lwd = 2)
  
  violins(data.frame(Gibbs.obj$Statistics[,3], MH.obj$Statistics[,3]), names = c("Gibbs", "M-H"), col = c("blue", "green"), connectcol = "transparent", main = "Cyclic Triads")
  abline(h = gergm.obj@stats[2,3], col = "red", lty = 2, lwd = 2)
  
  violins(data.frame(Gibbs.obj$Statistics[,4], MH.obj$Statistics[,4]), names = c("Gibbs", "M-H"), col = c("blue", "green"), connectcol = "transparent", main = "Reciprocity")
  abline(h = gergm.obj@stats[2,4], col = "red", lty = 2, lwd = 2)
  
  violins(data.frame(Gibbs.obj$Statistics[,5], MH.obj$Statistics[,5]), names = c("Gibbs", "M-H"), col = c("blue", "green"), connectcol = "transparent", main = "Transitive Triads")
  abline(h = gergm.obj@stats[2,5], col = "red", lty = 2, lwd = 2)
  
  violins(data.frame(Gibbs.obj$Statistics[,6], MH.obj$Statistics[,6]), names = c("Gibbs", "M-H"), col = c("blue", "green"), connectcol = "transparent", main = "Total Edge Weight")
  abline(h = gergm.obj@stats[2,6], col = "red", lty = 2, lwd = 2)
}
# Function to estimate gergms
gergm <- function(object, directed = c(TRUE, FALSE), MPLE.only = c(FALSE, TRUE),
                  transform.data = NULL,
                  method = c("Gibbs", "Metropolis"), 
                  max.num.iterations = 10, mc.num.iterations = 100, 
                  nsim = 500, thin = 1, 
                  shape.parameter = 1, weights = NULL, together = 0, 
                  MCMC.burnin = 100, seed = 123, tolerance = 0.01, gain.factor = 0) {
  ngpar <- dim(transform.data)[3] + 1
  set.seed(seed)
  alpha <- weights
  MPLE.only <- MPLE.only[1] #default is FALSE
  directed <- directed[1] #default is TRUE
  method <- method[1] #default is Gibbs
  res1 <- weighted.ergm.data(object, theta = NULL, alpha = alpha)
  statistics <- res1$statistics
  alphas <- res1$alphas
  net <- res1$net
  #cat("initial network for estimation")
  #print(head(net))
  possible.stats <- c("out2star", "in2star", "ctriads", "recip", "ttriads", 
                      "edgeweight")
  rhs.formula <- possible.stats[statistics > 0]
  rhs <- paste(rhs.formula, collapse = " + ")  #rewriting a formula for tnet
  
  # Flag if statistics do not meet requirement for Gibbs
  if (method == "Gibbs" & sum(alphas != 1) > 0) {
    warning(paste0("Some statistics do not have second order derivative = 0.",
                   " Switching to Metropolis"))
    method = "Metropolis"
  }
  
  # Flag if Metropolis is specified but Gibbs is OK
  if (method == "Metropolis" & sum(alphas != 1) == 0) {
    warning(paste0("All statistics have second order derivative = 0.",
                   " Consider switching to Gibbs for speed."))
    # method = "Gibbs"
  }
  
  # Estimation if a transformation is needed
  if (is.null(transform.data) != TRUE) {
    num.theta <- length(which(statistics > 0))
    gpar <- list()
    gpar$par <- c(mean(c(net)), rep(0, ngpar - 2), log(sd(c(net))))
    theta <- list()
    theta$par <- rep(0, num.theta)
    num.nodes <- nrow(net)
    
    # Alternately update lambda estimates and theta estimates
    for (i in 1:max.num.iterations) {
      ## If this is true, we don't need to simulate any networks, we simply give
      ## the MPLEs for the theta parameters at each iteration.
      if(MPLE.only == TRUE){
        cat("Updating Estimates -- Iteration:", i," \n")
        cat("Lambda Estimates", gpar$par,"\n")
        gpar.new <- optim(par = as.numeric(gpar$par), llg, formula = object, 
                          alpha = alpha, theta = as.numeric(theta$par), z = transform.data, 
                          method = "BFGS", together = together, 
                          hessian = T, control = list(fnscale = -1, trace = 6))
        cat("Lambda estimates", "\n")
        print(gpar.new$par)
        gpar.std.errors <- 1 / sqrt(abs(diag(gpar.new$hessian)))
        # Transform the unbounded weights to bounded weights via a t-distribution
        beta <- gpar.new$par[1:(length(gpar.new$par) - 1)]
        sig <- 0.01 + exp(gpar.new$par[length(gpar.new$par)])
        BZ <- 0
        for (j in 1:(dim(transform.data)[3])) {
          BZ <- BZ + beta[j] * transform.data[, , j]
        }
        net.new <<- pst(net, BZ, sig, 1) #make this a global variable (JW: not sure why we need this call..)
        
        # Rewrite the formula for net.new
        formula.new <- formula(paste0("net.new ~", rhs))
        
        # Estimate theta
        res1 <- weighted.ergm.data(formula.new, theta = theta$par, alpha = alpha)
        statistics <- res1$statistics
        net2 <- res1$net
        
        theta.new <- mple(net2, statistics = statistics, directed = directed)
        cat("theta.new", theta.new$par, "\n")
        cat("theta", theta$par, "\n")
        cat("statistics", statistics, "\n")
        theta.std.errors <- 1 / sqrt(abs(diag(theta.new$hessian)))
        
        if(i > 1){
          # Stop if lambda and theta estimates converge
          p.value1 <- rep(0,length(theta$par))
          count1 <- rep(0, length(theta$par))
          p.value2 <- rep(0, length(gpar$par))
          count2 <- rep(0, length(gpar$par))
          for(i in 1:length(theta$par)){
            #two sided z test
            p.value1[i] <- 2*pnorm(-abs((theta.new$par[i] - theta$par[i])/theta.std.errors[i]))
            #abs(theta.new$par[i] - theta$par[i]) > bounds[i]
            #if we reject any of the tests then convergence has not been reached!
            if(p.value1[i] < tolerance){count1[i] = 1}
          }
          for(i in 1:length(gpar$par)){
            #two sided z test
            p.value2[i] <- 2*pnorm(-abs((gpar.new$par[i] - gpar$par[i])/gpar.std.errors[i]))
            #if we reject any of the tests then convergence has not been reached!
            if(p.value2[i] < tolerance){count2[i] = 1}
          }
          cat("Theta p.values", "\n")
          print(p.value1)
          cat("Lambda p.values", "\n")
          print(p.value2)
          if (sum(count1) + sum(count2) == 0){
            message("Parameter estimates have converged")
            theta <- theta.new
            gpar <- gpar.new
            break
          }
        }
        theta <- theta.new
        gpar <- gpar.new
      }
      
      if(MPLE.only != TRUE){
        # Estimate lambda
        cat("Updating Estimates -- Iteration:", i," \n")
        cat("Lambda Estimates", gpar$par,"\n")
        gpar.new <- optim(par = as.numeric(gpar$par), llg, formula = object, 
                          alpha = alpha, theta = as.numeric(theta$par), z = transform.data, 
                          method = "BFGS", together = together, 
                          hessian = T, control = list(fnscale = -1, trace = 6))
        cat("Lambda estimates", "\n")
        print(gpar.new$par)
        gpar.std.errors <- 1 / sqrt(abs(diag(gpar.new$hessian)))
        # Transform the unbounded weights to bounded weights via a t-distribution
        beta <- gpar.new$par[1:(length(gpar.new$par) - 1)]
        sig <- 0.01 + exp(gpar.new$par[length(gpar.new$par)])
        BZ <- 0
        for (j in 1:(dim(transform.data)[3])) {
          BZ <- BZ + beta[j] * transform.data[, , j]
        }
        net.new <<- pst(net, BZ, sig, 1) #make this a global variable
        
        # Rewrite the formula for net.new
        formula.new <- formula(paste0("net.new ~", rhs))
        
        
        # Estimate theta
        theta.new <- mcmcmle(formula.obj = formula.new, 
                             num.draws = nsim,
                             mc.num.iterations = mc.num.iterations, 
                             thin = thin, MCMC.burnin = MCMC.burnin, 
                             theta = theta$par, alpha = alpha, 
                             directed = directed, method = method, 
                             shape.parameter = shape.parameter, together = together, 
                             tolerance = tolerance, seed2 = seed, gain.factor = gain.factor)
        
        # Calculate standard errors
        theta.std.errors <- 1 / sqrt(abs(diag(theta.new$hessian)))
        
        
        # Stop if lambda and theta estimates converge.
        # Convergence criterion is based on individual z-tests for each estimate
        p.value1 <- rep(0,length(theta$par))
        count1 <- rep(0, length(theta$par))
        p.value2 <- rep(0, length(gpar$par))
        count2 <- rep(0, length(gpar$par))
        for(i in 1:length(theta$par)){
          #two sided z test
          p.value1[i] <- 2*pnorm(-abs((theta.new$par[i] - theta$par[i])/theta.std.errors[i]))
          #abs(theta.new$par[i] - theta$par[i]) > bounds[i]
          #if we reject any of the tests then convergence has not been reached!
          if(p.value1[i] < tolerance){count1[i] = 1}
        }
        for(i in 1:length(gpar$par)){
          #two sided z test
          p.value2[i] <- 2*pnorm(-abs((gpar.new$par[i] - gpar$par[i])/gpar.std.errors[i]))
          #abs(theta.new$par[i] - theta$par[i]) > bounds[i]
          #if we reject any of the tests then convergence has not been reached!
          if(p.value2[i] < tolerance){count2[i] = 1}
        }
        cat("Theta p.values", "\n")
        print(p.value1)
        cat("Lambda p.values", "\n")
        print(p.value2)
        
        if (sum(count1) + sum(count2) == 0){
          message("Parameter estimates have converged")
          theta <- theta.new
          gpar <- gpar.new
          break
        }
        theta <- theta.new
        gpar <- gpar.new
      }
    }
    network <- net
    bounded.network <- net.new
    theta <- t(as.matrix(theta$par))
    theta <- rbind(theta, theta.std.errors)
    colnames(theta) <- rhs.formula
    rownames(theta) <- c("est", "se")
    theta <- as.data.frame(theta)
    lambda <- as.numeric(t(as.matrix(gpar$par)))
    lambda <- rbind(lambda, gpar.std.errors)
    lambda <- as.data.frame(lambda)
    rownames(lambda) <- c("est", "se")
    return(gergm.object(network, bounded.network, formula.obj, theta, lambda, 
                        alpha = alpha, together = together))
  }
  
  # Estimation if no transformation is needed
  if (is.null(transform.data) == TRUE) {
    num.theta <- length(which(statistics > 0))
    theta <- list()
    theta$par <- rep(0, num.theta)
    num.nodes <- nrow(net)
    if(MPLE.only == TRUE){
      res1 <- weighted.ergm.data(object, theta = theta$par, alpha = alpha)
      statistics <- res1$statistics
      net2 <- res1$net
      
      theta.new <- mple(net2, statistics = statistics, directed = directed)
      theta.std.errors <- 1 / sqrt(abs(diag(theta.new$hessian)))
      theta <- theta.new
      lambda <- as.data.frame(0)
    }
    
    if(MPLE.only != TRUE){
      theta.new <- mcmcmle(formula.obj = object, num.draws = nsim, 
                           mc.num.iterations = mc.num.iterations, thin = thin,
                           MCMC.burnin = MCMC.burnin, 
                           theta = theta$par, alpha = alpha, directed = directed, 
                           method = method, shape.parameter = shape.parameter, 
                           together = together, tolerance = tolerance, gain.factor = gain.factor)
      theta.std.errors <- 1 / sqrt(abs(diag(theta.new$hessian)))
      theta <- theta.new
      lambda <- as.data.frame(0)
    }
  }
  
  network <- net
  bounded.network <- net
  theta <- t(as.matrix(theta$par))
  theta <- rbind(theta, theta.std.errors)
  colnames(theta) <- rhs.formula
  rownames(theta) <- c("est", "se")
  theta <- as.data.frame(theta)
  return(gergm.object(network, bounded.network, formula.obj, theta, lambda, 
                      alpha = alpha, together = together))
}

