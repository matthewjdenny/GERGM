
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
                           possible.stats,
                           GERGM_Object,
                           force_x_theta_updates,
                           transformation_type) {

  # set the seed
  set.seed(seed)

  # set our exponential down weights
  alpha <- exponential_weights

  statistics <- GERGM_Object@stats_to_use
  alphas <- GERGM_Object@weights
  net <- GERGM_Object@network

  rhs.formula <- possible.stats[statistics > 0]
  rhs <- paste(rhs.formula, collapse = " + ")  #rewriting a formula for tnet

  # Flag if statistics do not meet requirement for Gibbs
  if (method == "Gibbs" & sum(GERGM_Object@weights != 1) > 0) {
    warning(paste0("Some statistics do not have second order derivative = 0.",
                   " Switching to Metropolis"))
    method = "Metropolis"
  }

  if(method == "Gibbs" & GERGM_Object@is_correlation_network){
    warning("Gibbs sampling is currently not implemented for correlation networks, switching to Metropolis Hastings.")
    method = "Metropolis"
  }

  # Flag if Metropolis is specified but Gibbs is OK
  if (method == "Metropolis" & sum(GERGM_Object@weights != 1) == 0) {
    warning(paste0("All statistics have second order derivative = 0.",
                   " Consider switching to Gibbs for speed."))
    # method = "Gibbs"
  }


  # Estimation if a transformation is needed
  if (length(GERGM_Object@data_transformation) > 0) {
    num.theta <- length(which(GERGM_Object@stats_to_use > 0))
    gpar <- list()
    gpar$par <- c(mean(c(GERGM_Object@network)),
                  rep(0, dim(GERGM_Object@data_transformation)[3] - 1),
                  log(sd(c(GERGM_Object@network))))
    theta <- list()
    theta$par <- rep(0, num.theta)
    num.nodes <- nrow(GERGM_Object@num_nodes)

    # Alternately update lambda estimates and theta estimates
    for (i in 1:max.num.iterations) {
      ## If this is true, we don't need to simulate any networks, we simply give
      ## the MPLEs for the theta parameters at each iteration.
      if(MPLE.only == TRUE){
        cat("Updating Estimates -- Iteration:", i," \n")
        GERGM_Object <- store_console_output(GERGM_Object,paste("Updating Estimates -- Iteration:", i," \n"))
        cat("Lambda Estimates", gpar$par,"\n")
        GERGM_Object <- store_console_output(GERGM_Object,paste("Lambda Estimates", gpar$par,"\n"))
        gpar.new <- optim(par = as.numeric(gpar$par),
                          llg,
                          alpha = GERGM_Object@weights,
                          theta = as.numeric(theta$par),
                          z = GERGM_Object@data_transformation,
                          method = "BFGS",
                          together = together,
                          GERGM_Object = GERGM_Object,
                          hessian = T,
                          control = list(fnscale = -1, trace = 6))
        cat("Lambda estimates", "\n")
        GERGM_Object <- store_console_output(GERGM_Object, "Lambda estimates\n")
        print(gpar.new$par)
        GERGM_Object <- store_console_output(GERGM_Object, toString(gpar.new$par))
        gpar.std.errors <- 1 / sqrt(abs(diag(gpar.new$hessian)))
        # Transform the unbounded weights to bounded weights via a t-distribution
        beta <- gpar.new$par[1:(length(gpar.new$par) - 1)]
        sig <- 0.01 + exp(gpar.new$par[length(gpar.new$par)])
        BZ <- 0
        for (j in 1:(dim(GERGM_Object@data_transformation)[3])) {
          BZ <- BZ + beta[j] * GERGM_Object@data_transformation[, , j]
        }

        if(transformation_type == "logcauchy" | transformation_type == "cauchy"){
          GERGM_Object@bounded.network <- pst(net, BZ, sig, 1)
        }
        if(transformation_type == "lognormal" | transformation_type == "gaussian"){
          GERGM_Object@bounded.network <- pst(net, BZ, sig, Inf)
        }

        num.nodes <- GERGM_Object@num_nodes
        triples <- t(combn(1:num.nodes, 3))
        pairs <- t(combn(1:num.nodes, 2))

        theta.new <- mple(GERGM_Object@bounded.network,
                          statistics = GERGM_Object@stats_to_use,
                          directed = directed)
        cat("theta.new", theta.new$par, "\n")
        GERGM_Object <- store_console_output(GERGM_Object,paste("theta.new", theta.new$par, "\n"))
        cat("theta", theta$par, "\n")
        GERGM_Object <- store_console_output(GERGM_Object,paste("theta", theta$par, "\n"))
        cat("statistics", statistics, "\n")
        GERGM_Object <- store_console_output(GERGM_Object,paste("statistics", statistics, "\n"))
        theta.std.errors <- 1 / sqrt(abs(diag(theta.new$hessian)))
        GERGM_Object@theta.par <- theta.new$par

        if(i > 1){
          # Stop if lambda and theta estimates converge
          p.value1 <- rep(0, length(theta$par))
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
          GERGM_Object <- store_console_output(GERGM_Object,paste("Theta p.values", "\n"))
          print(p.value1)
          GERGM_Object <- store_console_output(GERGM_Object,paste0(p.value1,collapse = " "))
          cat("Lambda p.values", "\n")
          GERGM_Object <- store_console_output(GERGM_Object,paste("Lambda p.values", "\n"))
          print(p.value2)
          GERGM_Object <- store_console_output(GERGM_Object,paste0(p.value2,collapse = " "))
          if (sum(count1) + sum(count2) == 0){
            message("Parameter estimates have converged")
            GERGM_Object <- store_console_output(GERGM_Object,"Parameter estimates have converged")
            GERGM_Object@theta_estimation_converged <- TRUE
            GERGM_Object@lambda_estimation_converged <- TRUE
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
        GERGM_Object <- store_console_output(GERGM_Object,paste("Updating Estimates -- Iteration:", i," \n"))
        cat("Lambda Estimates", gpar$par,"\n")
        GERGM_Object <- store_console_output(GERGM_Object,paste("Lambda Estimates", gpar$par,"\n"))
        gpar.new <- optim(par = as.numeric(gpar$par),
                          llg,
                          alpha = GERGM_Object@weights,
                          theta = as.numeric(theta$par),
                          z = GERGM_Object@data_transformation,
                          method = "BFGS",
                          together = together,
                          GERGM_Object = GERGM_Object,
                          hessian = T, control = list(fnscale = -1, trace = 6))
        cat("Lambda estimates", "\n")
        GERGM_Object <- store_console_output(GERGM_Object,"Lambda estimates\n")
        print(gpar.new$par)
        GERGM_Object <- store_console_output(GERGM_Object,toString(gpar.new$par))
        gpar.std.errors <- 1 / sqrt(abs(diag(gpar.new$hessian)))
        # Transform the unbounded weights to bounded weights via a t-distribution
        beta <- gpar.new$par[1:(length(gpar.new$par) - 1)]
        sig <- 0.01 + exp(gpar.new$par[length(gpar.new$par)])
        BZ <- 0
        for (j in 1:(dim(transform.data)[3])) {
          BZ <- BZ + beta[j] * transform.data[, , j]
        }

        if(transformation_type == "logcauchy" | transformation_type == "cauchy"){
          GERGM_Object@bounded.network <- pst(net, BZ, sig, 1)
        }
        if(transformation_type == "lognormal" | transformation_type == "gaussian"){
          GERGM_Object@bounded.network <- pst(net, BZ, sig, Inf)
        }

        num.nodes <- GERGM_Object@num_nodes
        triples <- t(combn(1:num.nodes, 3))
        pairs <- t(combn(1:num.nodes, 2))

        # Estimate theta
        ret_list <- MCMCMLE(num.draws = nsim,
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
							               possible.stats = possible.stats,
							               GERGM_Object= GERGM_Object,
							               force_x_theta_updates = force_x_theta_updates)

        theta.new <- ret_list[[1]]
        GERGM_Object <- ret_list[[2]]

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
          # break if we have hit model degeneracy
          if(GERGM_Object@theta_estimation_converged){
            if(p.value1[i] < tolerance){count1[i] = 1}
          }
        }
        for(i in 1:length(gpar$par)){
          #two sided z test
          p.value2[i] <- 2*pnorm(-abs((gpar.new$par[i] - gpar$par[i])/gpar.std.errors[i]))
          #abs(theta.new$par[i] - theta$par[i]) > bounds[i]
          #if we reject any of the tests then convergence has not been reached!
          if(p.value2[i] < tolerance){count2[i] = 1}
        }
        cat("Theta p.values", "\n")
        GERGM_Object <- store_console_output(GERGM_Object,paste("Theta p.values", "\n"))
        print(p.value1)
        GERGM_Object <- store_console_output(GERGM_Object,paste0(p.value1,collapse = " "))
        cat("Lambda p.values", "\n")
        GERGM_Object <- store_console_output(GERGM_Object,paste("Lambda p.values", "\n"))
        print(p.value2)
        GERGM_Object <- store_console_output(GERGM_Object,paste0(p.value2,collapse = " "))

        if (sum(count1) + sum(count2) == 0){
          message("Parameter estimates have converged")
          GERGM_Object <- store_console_output(GERGM_Object,"Parameter estimates have converged")
          GERGM_Object@theta_estimation_converged <- TRUE
          GERGM_Object@lambda_estimation_converged <- TRUE
          theta <- theta.new
          gpar <- gpar.new
          break
        }
        theta <- theta.new
        gpar <- gpar.new
      }
    }
    theta <- t(as.matrix(theta$par))
    theta <- rbind(theta, theta.std.errors)
    colnames(theta) <- rhs.formula
    rownames(theta) <- c("est", "se")
    theta <- as.data.frame(theta)
    lambda <- as.numeric(t(as.matrix(gpar$par)))
    lambda <- rbind(lambda, gpar.std.errors)
    lambda <- as.data.frame(lambda)
    rownames(lambda) <- c("est", "se")
    GERGM_Object@theta.coef <- theta
    GERGM_Object@lambda.coef <- lambda
    return(GERGM_Object)
  }

  # Estimation if no transformation is needed
  if (length(GERGM_Object@data_transformation) == 0) {
    GERGM_Object@lambda_estimation_converged <- TRUE
    num.theta <- length(which(GERGM_Object@stats_to_use > 0))
    theta <- list()
    theta$par <- rep(0, num.theta)
    num.nodes <- GERGM_Object@num_nodes
    if(MPLE.only == TRUE){
      ###
      statistics <- GERGM_Object@stats_to_use
      alphas <- GERGM_Object@weights
      cat("Estimating Theta via MPLE... \n")
      GERGM_Object <- store_console_output(GERGM_Object, "Estimating Theta via MPLE... \n")
      if(GERGM_Object@is_correlation_network){
        theta.init <- mple.corr(GERGM_Object@network, GERGM_Object@bounded.network,
                                statistics = GERGM_Object@stats_to_use,
                                directed = directed)
      }else{
        theta.init <- mple(GERGM_Object@bounded.network,
                           statistics = GERGM_Object@stats_to_use,
                           directed = directed)
      }

      #print(theta.init)
      GERGM_Object <- store_console_output(GERGM_Object,paste("\nMPLE Thetas: ", theta.init$par, "\n"))
      theta <- theta.init
      lambda <- as.data.frame(0)
      theta.std.errors <- 1 / sqrt(abs(diag(theta.init$hessian)))
      GERGM_Object@theta.par <- theta.init$par
    }

    if(MPLE.only != TRUE){
      ret_list <- MCMCMLE(num.draws = nsim,
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
                           possible.stats = possible.stats,
                           GERGM_Object = GERGM_Object,
                           force_x_theta_updates = force_x_theta_updates)

      theta.new <- ret_list[[1]]
      GERGM_Object <- ret_list[[2]]

      theta.std.errors <- 1 / sqrt(abs(diag(theta.new$hessian)))
      theta <- theta.new
      lambda <- as.data.frame(0)
    }
  }
  theta <- t(as.matrix(theta$par))
  theta <- rbind(theta, theta.std.errors)
  colnames(theta) <- rhs.formula
  rownames(theta) <- c("est", "se")
  theta <- as.data.frame(theta)
  GERGM_Object@theta.coef <- theta
  GERGM_Object@lambda.coef <- lambda
  return(GERGM_Object)
}

