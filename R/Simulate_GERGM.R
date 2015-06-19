# Simulate a gergm
Simulate_GERGM <- function(object, nsim, seed = NULL,
                           method , MCMC.burnin ,  
                           weights = object@weights,
                           coef = object@theta.coef, formula = object@formula,
                           thin , shape.parameter , together ,
                           seed1,
						   possible.stats ) {
  # object: an object of class "gergm"
  
  theta <- as.numeric(coef[1,])  # obtain the estimated thetas
  
  sample_every <- floor(1/thin)
  
  net <- object@bounded.network  # obtain the bounded network object
  alpha <- weights
  rhs <- names(object@theta.coef)[abs(theta) > 0]
  rhs <- paste(rhs, collapse = "+")
  formula.new <- as.formula(paste0("net ~ ", rhs))
  
  res1 <- Parse_Formula_Object(formula.new, possible.stats, theta = theta[abs(theta) > 0], 
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

    print("Thetas")
    print(thetas)
    
    samples <- Metropolis_Hastings_Sampler(
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
    
  }
  h.statistics = data.frame(out2stars = h.statistics[, 1], 
                            in2stars = h.statistics[, 2], 
                            ctriads = h.statistics[, 3], 
                            recip = h.statistics[, 4], 
                            ttriads = h.statistics[, 5], 
                            edgeweight = h.statistics[, 6])

  return(list(Networks = nets, Statistics = h.statistics,
              Acceptance.rate = acceptance.rate))
}