MCMCMLE <- function(num.draws,
                    mc.num.iterations,
                    tolerance,
                    thin = 1,
                    MCMC.burnin,
                    theta = NULL,
                    alpha = NULL,
                    directed ,
                    together ,
                    seed2 ,
                    gain.factor,
					          possible.stats,
					          GERGM_Object,
					          force_x_theta_updates,
					          verbose) {

  statistics <- GERGM_Object@stats_to_use
  alphas <- GERGM_Object@weights
  if(verbose){
    cat("Estimating Initial Values for Theta via MPLE... \n")
  }
  GERGM_Object <- store_console_output(GERGM_Object,"Estimating Initial Values for Theta via MPLE... \n")

  if(GERGM_Object@is_correlation_network){
    theta.init <- mple.corr(GERGM_Object@network, GERGM_Object@bounded.network,
                            statistics = GERGM_Object@stats_to_use,
                            directed = directed)
  }else{
    theta.init <- mple(GERGM_Object@bounded.network,
                       statistics = GERGM_Object@stats_to_use,
                       directed = directed)
  }
  if(verbose){
    cat("\nMPLE Thetas: ", theta.init$par, "\n")
  }
  GERGM_Object <- store_console_output(GERGM_Object, paste("\nMPLE Thetas: ", theta.init$par, "\n"))
  num.nodes <- GERGM_Object@num_nodes
  triples <- t(combn(1:num.nodes, 3))
  if(GERGM_Object@is_correlation_network){
    # initialize the network with the observed network
    initial_network <- GERGM_Object@network
    # calculate the statistics of the original network
    init.statistics <- h2(GERGM_Object@network,
                          triples = triples,
                          statistics = rep(1, length(possible.stats)),
                          alphas = alphas, together = together)
    obs.stats <- h2(GERGM_Object@network,
                    triples = triples,
                    statistics = GERGM_Object@stats_to_use,
                    alphas = alphas,
                    together = together)
  }else{
    # initialize the network with the observed network
    initial_network <- GERGM_Object@bounded.network
    # calculate the statistics of the original network
    init.statistics <- h2(GERGM_Object@bounded.network,
                          triples = triples,
                          statistics = rep(1, length(possible.stats)),
                          alphas = alphas, together = together)
    obs.stats <- h2(GERGM_Object@bounded.network,
                    triples = triples,
                    statistics = GERGM_Object@stats_to_use,
                    alphas = alphas,
                    together = together)
  }


  #cat("Observed Values of Selected Statistics:", "\n", obs.stats, "\n")
  ####################################################################
  alps <- alphas[which(statistics == 1)]
  GERGM_Object@reduced_weights <- alps
  GERGM_Object@theta.par <- theta.init$par

  # if we are not doing a fisher update
  theta <- list()
  theta$par <- theta.init$par

  # if we are going to do a fisher update to MPLE thetas
  if(gain.factor > 0){
    GERGM_Object <- Simulate_GERGM(GERGM_Object,
                                   nsim = ceiling(20/thin),
                                   together = together,
                                   thin = thin,
                                   MCMC.burnin = MCMC.burnin,
                                   seed1 = seed2,
                                   possible.stats = possible.stats,
                                   verbose = verbose)

    hsn <- GERGM_Object@MCMC_output$Statistics[,which(GERGM_Object@stats_to_use == 1)]

    #Calculate covariance estimate (to scale initial guess theta.init)
    z.bar <- NULL
    if(class(hsn) == "numeric"){
      hsn <- matrix(hsn,ncol =1,nrow = length(hsn))
      z.bar <- sum(hsn) / 20
    }else{
      z.bar <- colSums(hsn) / 20
    }

    #cat("z.bar", "\n", z.bar, "\n")
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
    if(verbose){
      cat("Adjusted Initial Thetas After Fisher Update:",theta$par, "\n\n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,paste("Adjusted Initial Thetas After Fisher Update:",theta$par, "\n\n"))
  }

  ##########################################################################
  ## Simulate new networks
  for (i in 1:mc.num.iterations) {
    GERGM_Object@theta.par <- as.numeric(theta$par)
    GERGM_Object <- Simulate_GERGM(GERGM_Object,
                           nsim = num.draws,
                           together = together,
                           thin = thin,
                           MCMC.burnin = MCMC.burnin,
                           seed1 = seed2,
                           possible.stats = possible.stats,
                           verbose = verbose)

    hsn <- GERGM_Object@MCMC_output$Statistics[,which(statistics == 1)]
    hsn.tot <- GERGM_Object@MCMC_output$Statistics

    # deal with case where we only have one statistic
    if(class(hsn.tot) == "numeric"){
      hsn.tot <- matrix(hsn.tot,ncol =1,nrow = length(hsn.tot))
      stats.data <- data.frame(Observed = init.statistics,
                               Simulated = mean(hsn.tot))
    }else{
      stats.data <- data.frame(Observed = init.statistics,
                               Simulated = colMeans(hsn.tot))
    }

    rownames(stats.data) <- possible.stats
    cat("Simulated (averages) and observed network statistics...\n")
    print(stats.data)
    GERGM_Object <- store_console_output(GERGM_Object,toString(stats.data))
    if(verbose){
      cat("\nOptimizing theta estimates... \n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,"\nOptimizing Theta Estimates... \n")
    if(verbose){
    theta.new <- optim(par = theta$par,
                       log.l,
                       alpha = GERGM_Object@reduced_weights,
                       hsnet = hsn,
                       ltheta = as.numeric(theta$par),
                       together = together,
                       possible.stats= possible.stats,
                       GERGM_Object = GERGM_Object,
                       method = "BFGS",
                       hessian = T,
                       control = list(fnscale = -1, trace = 6))
    }else{
      theta.new <- optim(par = theta$par,
                         log.l,
                         alpha = GERGM_Object@reduced_weights,
                         hsnet = hsn,
                         ltheta = as.numeric(theta$par),
                         together = together,
                         possible.stats= possible.stats,
                         GERGM_Object = GERGM_Object,
                         method = "BFGS",
                         hessian = T,
                         control = list(fnscale = -1, trace = 0))
    }
    if(verbose){
      cat("\n", "Theta Estimates: ", paste0(theta.new$par,collapse = " "), "\n",sep = "")
    }
    GERGM_Object <- store_console_output(GERGM_Object,paste("\n", "Theta Estimates: ", paste0(theta.new$par,collapse = " "), "\n",sep = ""))
    theta.std.errors <- 1 / sqrt(abs(diag(theta.new$hessian)))
    # Calculate the p-value based on a z-test of differences
    # The tolerance is the alpha at which differences are significant
    p.value <- rep(0,length(as.numeric(theta$par)))
    count <- rep(0, length(as.numeric(theta$par)))
    for(j in 1:length(theta$par)){
      #two sided z test
      p.value[j] <- 2*pnorm(-abs((as.numeric(theta.new$par)[j] - as.numeric(theta$par)[j])/theta.std.errors[j]))
      #abs(theta.new$par[i] - theta$par[i]) > bounds[i]
      #if we reject any of the tests then convergence has not been reached!
      if(p.value[j] < tolerance){count[j] = 1}
    }
    if(verbose){
      cat("\np.values for two-sided z-test of difference between current and updated theta estimates:\n\n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,"\np.values for two-sided z-test of difference between current and updated theta estimates:\n\n")
    if(verbose){
      cat(round(p.value,3), "\n \n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,paste(p.value, "\n \n"))

    if(max(abs(theta.new$par)) > 10000000){
      message("Parameter estimates appear to have become degenerate, returning previous thetas. Model output should not be trusted. Try specifying a larger number of simulations or a different parameterization.")
      GERGM_Object <- store_console_output(GERGM_Object,"Parameter estimates appear to have become degenerate, returning previous thetas. Model output should not be trusted. Try specifying a larger number of simulations or a different parameterization.")
      return(list(theta.new,GERGM_Object))
    }

    # check to see if we had a zero percent accept rate if using MH, and if so,
    # then adjust proposal variance and try again -- do not signal convergence.
    allow_convergence = TRUE
    if (GERGM_Object@estimation_method == "Metropolis") {
      if (GERGM_Object@MCMC_output$Acceptance.rate == 0){
        cat("Acceptance rate was zero, decreasing proposal variance...\n")
        allow_convergence = FALSE
      }
    }
    if (allow_convergence) {
      if (sum(count) == 0){
        #conditional to check and see if we are requiring more updates
        if(i >= force_x_theta_updates){
          if(verbose){
            message("Parameter estimates have converged")
          }
          GERGM_Object <- store_console_output(GERGM_Object,"Parameter estimates have converged")
          GERGM_Object@theta_estimation_converged <- TRUE
          return(list(theta.new,GERGM_Object))
        }else{
          if(verbose){
            message(paste("Forcing",force_x_theta_updates,"iterations of theta updates..."),sep = " ")
          }
          GERGM_Object <- store_console_output(GERGM_Object,paste("Forcing",force_x_theta_updates,"iterations of theta updates..."))
        }
      }
    }
    theta <- theta.new
    GERGM_Object@theta.par <- as.numeric(theta$par)
  }
  return(list(theta.new,GERGM_Object))
}
