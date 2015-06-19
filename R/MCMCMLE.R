MCMCMLE <- function(formula.obj,
                    num.draws,
                    mc.num.iterations,
                    tolerance,
                    thin = 1,
                    MCMC.burnin,
                    theta = NULL,
                    alpha = NULL,
                    directed = c(TRUE, FALSE),
                    method = c("Gibbs", "Metropolis"),
                    shape.parameter = 1,
                    take.sample.every = 1,
                    together = 0,
                    seed2 = 10000,
                    gain.factor) {
  directed <- directed[1]
  method <- method[1]
  res1 <- weighted.ergm.data(formula.obj, theta = theta, alpha = alpha)
  statistics <- res1$statistics
  alphas <- res1$alphas
  #cat("alphas in mcmcmle")
  #print(alphas)
  net2 <- res1$net

  theta.init <- mple(net2, statistics = statistics, directed = directed)

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
