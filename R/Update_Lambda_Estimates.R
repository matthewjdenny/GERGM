Update_Lambda_Estimates <- function(i,
                                    gpar,
                                    theta,
                                    together,
                                    verbose,
                                    net,
                                    GERGM_Object){
  # if we are usinga correlation network, do beta regression
  if (GERGM_Object@is_correlation_network) {
    stop("Currently not implemented! Set omit_intercept_term = TRUE and include an edges term in specification...")
  } else if (GERGM_Object@beta_correlation_model) {
    ### This si where the code from James starts -- this is a conditional if we
    ### are using the beta model

    # initialize everything
    temp <- net[lower.tri(net)] #this will be n choose 2 in length
    beta <- rep(0, dim(GERGM_Object@data_transformation)[3])
    phi <- 1
    gpar$par <- c(beta, phi)

    # give the user some information
    cat("Updating Estimates -- Iteration:", i," \n")
    GERGM_Object <- store_console_output(GERGM_Object,paste("Updating Estimates -- Iteration:", i," \n"))
    if (verbose) {
      cat("Lambda Estimates (from Beta model)", gpar$par,"\n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,
      paste("Lambda Estimates (from Beta model)", gpar$par,"\n"))
    gpar.new <- NULL

    # run optime to update the parameters from the beta regression
    if (verbose) {
      gpar.new <- optim(par = as.numeric(gpar$par),
                        llg.beta,
                        theta = as.numeric(theta$par),
                        method = "BFGS",
                        together = together,
                        GERGM_Object = GERGM_Object,
                        hessian = T,
                        control = list(fnscale = -1, trace = 6))
    } else {
      gpar.new <- optim(par = as.numeric(gpar$par),
                        llg.beta,
                        theta = as.numeric(theta$par),
                        method = "BFGS",
                        together = together,
                        GERGM_Object = GERGM_Object,
                        hessian = T,
                        control = list(fnscale = -1, trace = 0))
    }

    #Transform to bounded network X via beta cdf
    beta <- gpar.new$par[1:(length(gpar.new$par) - 1)]
    mu <- logistic(beta, GERGM_Object@data_transformation)
    phi <- gpar.new$par[length(gpar.new$par)]
    shape1 <- mu * phi
    shape2 <- (1 - mu) * phi
    X <- pbeta(lower.tri((GERGM_Object@bounded.network + 1)/2),
               shape1 = shape1, shape2 = shape2)

    # store our new bounded network
    GERGM_Object@bounded.network <- X

    # tell the user abotu our new estimates
    if (verbose) {
      cat("Lambda estimates (from Beta model) \n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,
                      "Lambda estimates (from Beta model)\n")
    if (verbose) {
      print(gpar.new$par)
    }

    # get standard errors
    GERGM_Object <- store_console_output(GERGM_Object, toString(gpar.new$par))
    temp <- calculate_standard_errors(hessian = gpar.new$hessian,
                                      GERGM_Object = GERGM_Object)
    gpar.std.errors <- temp$std_errors
    GERGM_Object <- temp$GERGM_Object

    # store mu and phi for later use in reverse transformation
    GERGM_Object@BZ <- mu
    GERGM_Object@BZstdev <- phi

  }  else {
    #do our normal t regression
    cat("Updating Estimates -- Iteration:", i," \n")
    GERGM_Object <- store_console_output(GERGM_Object,paste("Updating Estimates -- Iteration:", i," \n"))
    if (verbose) {
      cat("Lambda Estimates", gpar$par,"\n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,paste("Lambda Estimates", gpar$par,"\n"))
    gpar.new <- NULL
    if (verbose) {
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
    }else{
      gpar.new <- optim(par = as.numeric(gpar$par),
                        llg,
                        alpha = GERGM_Object@weights,
                        theta = as.numeric(theta$par),
                        z = GERGM_Object@data_transformation,
                        method = "BFGS",
                        together = together,
                        GERGM_Object = GERGM_Object,
                        hessian = T,
                        control = list(fnscale = -1, trace = 0))
    }
    if (verbose) {
      cat("Lambda estimates", "\n")
    }
    GERGM_Object <- store_console_output(GERGM_Object, "Lambda estimates\n")
    if (verbose) {
      print(gpar.new$par)
    }
    GERGM_Object <- store_console_output(GERGM_Object, toString(gpar.new$par))
    temp <- calculate_standard_errors(hessian = gpar.new$hessian,
                                      GERGM_Object = GERGM_Object)
    gpar.std.errors <- temp$std_errors
    GERGM_Object <- temp$GERGM_Object
    # Transform the unbounded weights to bounded weights via a t-distribution
    beta <- gpar.new$par[1:(length(gpar.new$par) - 1)]
    sig <- 0.01 + exp(gpar.new$par[length(gpar.new$par)])
    BZ <- 0
    for (j in 1:(dim(GERGM_Object@data_transformation)[3])) {
      BZ <- BZ + beta[j] * GERGM_Object@data_transformation[, , j]
    }

    #store so we can transform back
    GERGM_Object@BZ <- BZ
    GERGM_Object@BZstdev <- sig
    transformation_type <- GERGM_Object@transformation_type
    if (transformation_type == "logcauchy") {
      GERGM_Object@bounded.network <- pst(log(net), BZ, sig, 1)
    }
    if (transformation_type == "cauchy") {
      GERGM_Object@bounded.network <- pst(net, BZ, sig, 1)
    }
    if (transformation_type == "lognormal") {
      GERGM_Object@bounded.network <- pst(log(net), BZ, sig, Inf)
    }
    if (transformation_type == "gaussian") {
      GERGM_Object@bounded.network <- pst(net, BZ, sig, Inf)
    }
  } # end of conditional if we are doing standard t regression
  return(list(GERGM_Object = GERGM_Object,
              gpar.new = gpar.new,
              gpar.std.errors  = gpar.std.errors))
}


