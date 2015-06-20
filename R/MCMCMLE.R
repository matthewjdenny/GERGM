MCMCMLE <- function(formula.obj,
                    num.draws,
                    mc.num.iterations,
                    tolerance,
                    thin = 1,
                    MCMC.burnin,
                    theta = NULL,
                    alpha = NULL,
                    directed ,
                    method ,
                    shape.parameter ,
                    take.sample.every ,
                    together ,
                    seed2 ,
                    gain.factor,
					          possible.stats) {

  res1 <- Parse_Formula_Object(formula.obj, possible.stats, theta = theta, alpha = alpha)
  statistics <- res1$statistics
  alphas <- res1$alphas
  net2 <- res1$net

  theta.init <- mple(net2, statistics = statistics, directed = directed)

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

  object <- Create_GERGM_Object_From_Formula(formula.obj,
                                             theta.coef = theta.init$par,
                                             possible.stats = possible.stats,
                                             weights = alps,
                                             together = together)
  temp <- Simulate_GERGM(object,
                         nsim = ceiling(20/thin),
                         method = method,
                         shape.parameter = shape.parameter,
                         together = together,
                         thin = thin,
                         MCMC.burnin = MCMC.burnin,
                         seed1 = seed2,
                         possible.stats = possible.stats)

  hsn <- temp$Statistics[,which(statistics == 1)]

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
    object <- Create_GERGM_Object_From_Formula(formula.obj,
                                               theta.coef = theta$par,
                                               possible.stats = possible.stats,
                                               weights = alps,
                                               together = together)
    temp <- Simulate_GERGM(object,
                           nsim = num.draws,
                           method = method,
                           shape.parameter = shape.parameter,
                           together = together,
                           thin = thin,
                           MCMC.burnin = MCMC.burnin,
                           seed1 = seed2,
                           possible.stats = possible.stats)


    #just use what gets returned
    hsn <- temp$Statistics[,which(statistics == 1)]


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

    theta.new <- optim(par = theta$par,
                       log.l,
                       alpha = alps,
                       formula = formula.obj,
                       hsnet = hsn,
                       ltheta = as.numeric(theta$par),
                       together = together,
                       possible.stats= possible.stats,
                       method = "BFGS",
                       hessian = T,
                       control = list(fnscale = -1, trace = 5))
    print(paste0("Theta = ", theta.new$par))
    theta.std.errors <- 1 / sqrt(abs(diag(theta.new$hessian)))
    # Calculate the p-value based on a z-test of differences
    # The tolerance is the alpha at which differences are significant
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
