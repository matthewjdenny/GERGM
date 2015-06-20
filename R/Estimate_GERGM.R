
# Function to estimate gergms
Estimate_GERGM <- function(formula_object,
                           directed,
                           MPLE.only,
                           transform.data = NULL,
                           method,
                           max.num.iterations,
                           mc.num.iterations,
                           nsim,
                           thin,
                           shape.parameter,
                           exponential_weights = NULL,
                           together,
                           MCMC.burnin,
                           seed,
                           tolerance,
                           gain.factor,
                           possible.stats) {

  #' set the seed
  set.seed(seed)

  #' set our exponential down weights
  alpha <- exponential_weights

  coef(Model_Output)

  #' parse the formula object into the network, its raw statistics,
  #' alpha weights, and theta parameters
  Parsed_Formula_Object <- Parse_Formula_Object(formula_object,
                                                possible.stats,
                                                theta = NULL,
                                                alpha = alpha)

  #' extract raw network statistics, alphas and the sociomatrix from the
  #' Parsed_Formula_Object
  statistics <- Parsed_Formula_Object$statistics
  alphas <- Parsed_Formula_Object$alphas
  net <- Parsed_Formula_Object$net

  #cat("initial network for estimation")
  #print(head(net))
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
    gpar$par <- c(mean(c(net)), rep(0, dim(transform.data)[3] -3), log(sd(c(net))))
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
        gpar.new <- optim(par = as.numeric(gpar$par), llg, formula = formula_object,
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
        res1 <- Parse_Formula_Object(formula.new,
                                     possible.stats,
                                     theta = theta$par,
                                     alpha = alpha)
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
        gpar.new <- optim(par = as.numeric(gpar$par), llg, formula = formula_object,
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
        theta.new <- MCMCMLE(formula.obj = formula.new,
                             num.draws = nsim,
                             mc.num.iterations = mc.num.iterations,
                             thin = thin, MCMC.burnin = MCMC.burnin,
                             theta = theta$par,
                             alpha = alpha,
                             directed = directed,
                             method = method,
                             shape.parameter = shape.parameter,
                             together = together,
                             tolerance = tolerance,
                             seed2 = seed,
                             gain.factor = gain.factor,
							               possible.stats = possible.stats)

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
    return(Create_GERGM_Object(network,
                               bounded.network,
                               formula.obj,
                               theta,
                               lambda,
                               alpha = alpha,
                               together = together,
                               possible.stats))
  }

  # Estimation if no transformation is needed
  if (is.null(transform.data) == TRUE) {
    num.theta <- length(which(statistics > 0))
    theta <- list()
    theta$par <- rep(0, num.theta)
    num.nodes <- nrow(net)
    if(MPLE.only == TRUE){
      res1 <- Parse_Formula_Object(formula_object,
                                   possible.stats,
                                   theta = theta$par,
                                   alpha = alpha)
      statistics <- res1$statistics
      net2 <- res1$net

      theta.new <- mple(net2, statistics = statistics, directed = directed)
      theta.std.errors <- 1 / sqrt(abs(diag(theta.new$hessian)))
      theta <- theta.new
      lambda <- as.data.frame(0)
    }

    if(MPLE.only != TRUE){
      theta.new <- MCMCMLE(formula.obj = formula_object,
                           num.draws = nsim,
                           mc.num.iterations = mc.num.iterations,
                           thin = thin,
                           MCMC.burnin = MCMC.burnin,
                           theta = theta$par,
                           alpha = alpha,
                           directed = directed,
                           method = method,
                           shape.parameter = shape.parameter,
                           together = together,
                           tolerance = tolerance,
                           seed2 = seed,
                           gain.factor = gain.factor,
                           possible.stats = possible.stats)

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
  Model_Output@theta.coef <<- theta
  Model_Output@lambda.coef <<- lambda
#   return(Create_GERGM_Object(network,
#                              bounded.network,
#                              formula_object,
#                              theta,
#                              lambda,
#                              alpha = alpha,
#                              together = together,
#                              possible.stats))
}

