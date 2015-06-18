
# Function to estimate gergms
Estimate_GERGM <- function(object, directed = c(TRUE, FALSE), MPLE.only = c(FALSE, TRUE),
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

