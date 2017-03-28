
# Function to estimate gergms
Estimate_GERGM <- function(formula_object,
                           MPLE.only,
                           max.num.iterations,
                           mc.num.iterations,
                           seed,
                           tolerance,
                           possible.stats,
                           GERGM_Object,
                           force_x_theta_updates,
                           transformation_type,
                           verbose,
                           force_x_lambda_updates,
                           stop_for_degeneracy) {

  # set the seed
  set.seed(seed)

  # Flag if statistics do not meet requirement for Gibbs
  if (GERGM_Object@estimation_method == "Gibbs") {
    warning("Gibbs is currently not supported, considder switching to MH...")
  }

  # Estimation if a transformation is needed
  if (length(GERGM_Object@data_transformation) > 0) {
    num.theta <- length(GERGM_Object@stats_to_use)
    gpar <- list()

    # initialize gpar for lambda/beta optimization
    if (GERGM_Object@beta_correlation_model) {
      beta <- rep(0.5, dim(GERGM_Object@data_transformation)[3])
      phi <- 1
      gpar$par <- c(beta, phi)
    } else {
      gpar$par <- c(rep(0, dim(GERGM_Object@data_transformation)[3]),
                    log(sd(c(GERGM_Object@network))))
      # find the intercept term and set it equal to the mean of the network
      intercept <- mean(c(GERGM_Object@network))
      index <- which(dimnames(GERGM_Object@data_transformation)[[3]] == "intercept")
      if (length(index) == 1) {
        gpar$par[index] <- intercept
      }
      # give names to the lambda estimates
      names(gpar$par) <- c(dimnames(GERGM_Object@data_transformation)[[3]],"dispersion parameter")
    }
    theta <- list()
    theta$par <- rep(0, num.theta)
    num.nodes <- nrow(GERGM_Object@num_nodes)

    # Alternately update lambda estimates and theta estimates
    for (i in 1:max.num.iterations) {
      ## If this is true, we don't need to simulate any networks, we simply give
      ## the MPLEs for the theta parameters at each iteration.
      if(MPLE.only == TRUE){
        updates <- Update_Lambda_Estimates(
          i = i,
          gpar = gpar,
          theta = theta,
          together = GERGM_Object@downweight_statistics_together,
          verbose = verbose,
          net = GERGM_Object@network,
          GERGM_Object = GERGM_Object)

        GERGM_Object <- updates$GERGM_Object
        gpar.new <- updates$gpar.new
        gpar.std.errors <- updates$gpar.std.errors

        # get MPLE thetas
        MPLE_Results <- run_mple(GERGM_Object = GERGM_Object,
                                 verbose = verbose,
                                 seed2 = seed,
                                 possible.stats = possible.stats)

        GERGM_Object <- MPLE_Results$GERGM_Object
        theta.new <- MPLE_Results$theta
        statistics <- MPLE_Results$statistics
        init.statistics <- MPLE_Results$init.statistics

        if (GERGM_Object@covariate_terms_only) {
          cat("Forcing edges parameter to zero since it was not specified...\n")
          theta.new$par <- 0
          statistics <- init.statistics
        }

        if (verbose) {
          cat("theta.new", theta.new$par, "\n")
        }
        GERGM_Object <- store_console_output(GERGM_Object,
          paste("theta.new", theta.new$par, "\n"))
        if (verbose) {
          cat("theta", theta$par, "\n")
        }
        GERGM_Object <- store_console_output(GERGM_Object,
          paste("theta", theta$par, "\n"))
        if (verbose) {
          cat("statistics", GERGM_Object@stats_to_use, "\n")
        }
        GERGM_Object <- store_console_output(GERGM_Object,
          paste("statistics", GERGM_Object@stats_to_use, "\n"))

        temp <- calculate_standard_errors(hessian = MPLE_Results$hessian,
                                          GERGM_Object = GERGM_Object)
        theta.std.errors <- temp$std_errors
        GERGM_Object <- temp$GERGM_Object
        # theta.std.errors <- 1 / sqrt(abs(diag(theta.new$hessian)))
        print(theta.new$par)
        GERGM_Object@theta.par <- theta.new$par

        if (GERGM_Object@covariate_terms_only) {
          theta.std.errors <- 0.00001
        }

        if (i > 1) {
          # Stop if lambda and theta estimates converge
          p.value1 <- rep(0, length(theta$par))
          count1 <- rep(0, length(theta$par))
          p.value2 <- rep(0, length(gpar$par))
          count2 <- rep(0, length(gpar$par))
          for (i in 1:length(theta$par)) {
            #two sided z test
            p.value1[i] <- 2 * pnorm(-abs((theta.new$par[i] -
                                           theta$par[i]) / theta.std.errors[i]))
            #abs(theta.new$par[i] - theta$par[i]) > bounds[i]
            #if we reject any of the tests then convergence has not been reached!
            if (p.value1[i] < tolerance) {count1[i] = 1}
          }
          for (i in 1:length(gpar$par)) {
            #two sided z test
            p.value2[i] <- 2 * pnorm(-abs((gpar.new$par[i] -
                                           gpar$par[i]) / gpar.std.errors[i]))
            #if we reject any of the tests then convergence has not been reached!
            if (p.value2[i] < tolerance) {count2[i] = 1}
          }
          if (verbose) {
            cat("Theta p-values", "\n")
          }
          GERGM_Object <- store_console_output(GERGM_Object,
                                               paste("Theta p.values", "\n"))
          if (verbose) {
            names(p.value1) <- colnames(GERGM_Object@theta.coef)
            print(p.value1)
          }
          GERGM_Object <- store_console_output(GERGM_Object,
                                               paste0(p.value1,collapse = " "))
          if (verbose) {
            cat("Lambda p-values", "\n")
          }
          GERGM_Object <- store_console_output(GERGM_Object,
                                               paste("Lambda p.values", "\n"))
          if (verbose) {
            names(p.value2) <- c(dimnames(GERGM_Object@data_transformation)[[3]],"dispersion parameter")
            print(p.value2)
          }
          GERGM_Object <- store_console_output(GERGM_Object,
                                               paste0(p.value2,collapse = " "))
          if (sum(count1) + sum(count2) == 0) {
            message("Theta parameter estimates have converged...")
            if (force_x_lambda_updates > i) {
              cat("Forcing",force_x_lambda_updates,"lambda updates...\n")
            } else {
              GERGM_Object <- store_console_output(GERGM_Object,
                "Parameter estimates have converged")
              GERGM_Object@theta_estimation_converged <- TRUE
              GERGM_Object@lambda_estimation_converged <- TRUE
              theta <- theta.new
              gpar <- gpar.new
              break
            }
          }
        }
        theta <- theta.new
        gpar <- gpar.new
      }

      if (MPLE.only != TRUE) {
        # Estimate lambda
        updates <- Update_Lambda_Estimates(
          i = i,
          gpar = gpar,
          theta = theta,
          together = GERGM_Object@downweight_statistics_together,
          verbose = verbose,
          net = GERGM_Object@network,
          GERGM_Object = GERGM_Object)

        GERGM_Object <- updates$GERGM_Object
        gpar.new <- updates$gpar.new
        gpar.std.errors <- updates$gpar.std.errors

        if (GERGM_Object@covariate_terms_only) {
          theta.new <- list()
          theta.new$par <- 0
          theta.std.errors <- 0.00001
        } else {
          # Estimate theta
          ret_list <- MCMCMLE(
            mc.num.iterations = mc.num.iterations,
            theta = theta$par,
            tolerance = tolerance,
            seed2 = seed,
            possible.stats = possible.stats,
            GERGM_Object = GERGM_Object,
            force_x_theta_updates = force_x_theta_updates,
            verbose = verbose,
            outter_iteration_number = i,
            stop_for_degeneracy = stop_for_degeneracy)

          if (GERGM_Object@using_slackr_integration) {
            time <- Sys.time()
            message <- paste("Theta parameter estimates converged at:",
                             toString(time))
            slackr::slackr_bot(
              message,
              channel = GERGM_Object@slackr_integration_list$channel,
              username = GERGM_Object@slackr_integration_list$model_name,
              incoming_webhook_url = GERGM_Object@slackr_integration_list$incoming_webhook_url)
          }

          theta.new <- ret_list[[1]]
          GERGM_Object <- ret_list[[2]]

          # Calculate standard errors
          temp <- calculate_standard_errors(hessian = theta.new$hessian,
                                            GERGM_Object = GERGM_Object)
          theta.std.errors <- temp$std_errors
          GERGM_Object <- temp$GERGM_Object
        }

        # Stop if lambda and theta estimates converge.
        # Convergence criterion is based on individual z-tests for each estimate
        p.value2 <- rep(0, length(gpar$par))
        count2 <- rep(0, length(gpar$par))

        for (i in 1:length(gpar$par)) {
          #two sided z test
          p.value2[i] <- 2*pnorm(-abs((gpar.new$par[i] - gpar$par[i]) /
                                        gpar.std.errors[i]))
          #if we reject any of the tests then convergence has not been reached!
          if (p.value2[i] < tolerance) {count2[i] = 1}
        }

        if (verbose) {
          cat("Lambda p.values", "\n")
        }
        GERGM_Object <- store_console_output(GERGM_Object,
                                             paste("Lambda p.values", "\n"))
        if (verbose) {
          names(p.value2) <- c(dimnames(GERGM_Object@data_transformation)[[3]],"dispersion parameter")
          print(p.value2)
        }
        GERGM_Object <- store_console_output(GERGM_Object,
                                             paste0(p.value2,collapse = " "))

        if (GERGM_Object@using_slackr_integration) {
          time <- Sys.time()
          message <- paste("Lambda parameter estimate p-values at:",
                           toString(time))
          p_values <- paste(round(p.value2,3))
          slackr::slackr_bot(
            message,
            p_values,
            channel = GERGM_Object@slackr_integration_list$channel,
            username = GERGM_Object@slackr_integration_list$model_name,
            incoming_webhook_url = GERGM_Object@slackr_integration_list$incoming_webhook_url)
        }

        if (sum(count2) == 0) {
          message("Lambda parameter estimates have converged...")
          if (force_x_lambda_updates > i) {
            message("Forcing",force_x_lambda_updates,"lambda updates...\n")
          } else {
            GERGM_Object <- store_console_output(GERGM_Object,
                              "Parameter estimates have converged")
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
    }
    theta <- t(as.matrix(theta$par))
    theta <- rbind(theta, theta.std.errors)
    colnames(theta) <- GERGM_Object@theta_names
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
    if (MPLE.only == TRUE) {
      if (verbose) {
        cat("Estimating Theta via MPLE... \n")
      }
      GERGM_Object <- store_console_output(GERGM_Object,
                                           "Estimating Theta via MPLE... \n")
      # get MPLE thetas
      MPLE_Results <- run_mple(GERGM_Object = GERGM_Object,
                               verbose = verbose,
                               seed2 = seed,
                               possible.stats = possible.stats)

      GERGM_Object <- MPLE_Results$GERGM_Object
      theta.init <- MPLE_Results$theta
      statistics <- MPLE_Results$statistics
      init.statistics <- MPLE_Results$init.statistics

      GERGM_Object <- store_console_output(GERGM_Object,paste("\nMPLE Thetas: ",
                                                              theta.init$par,
                                                              "\n"))
      theta <- theta.init
      lambda <- as.data.frame(0)

      temp <- calculate_standard_errors(hessian = MPLE_Results$hessian,
                                        GERGM_Object = GERGM_Object)
      theta.std.errors <- temp$std_errors
      GERGM_Object <- temp$GERGM_Object
      GERGM_Object@theta.par <- theta.init$par
    }

    if (MPLE.only != TRUE) {
      ret_list <- MCMCMLE(
        mc.num.iterations = mc.num.iterations,
        theta = theta$par,
        tolerance = tolerance,
        seed2 = seed,
        possible.stats = possible.stats,
        GERGM_Object = GERGM_Object,
        force_x_theta_updates = force_x_theta_updates,
        verbose = verbose,
        stop_for_degeneracy = stop_for_degeneracy)

      theta.new <- ret_list[[1]]
      GERGM_Object <- ret_list[[2]]

      temp <- calculate_standard_errors(hessian = theta.new$hessian,
                                        GERGM_Object = GERGM_Object)
      theta.std.errors <- temp$std_errors
      GERGM_Object <- temp$GERGM_Object
      theta <- theta.new
      lambda <- as.data.frame(0)
    }
  }
  theta <- t(as.matrix(theta$par))
  theta <- rbind(theta, theta.std.errors)
  colnames(theta) <- GERGM_Object@theta_names
  rownames(theta) <- c("est", "se")
  theta <- as.data.frame(theta)
  GERGM_Object@theta.coef <- theta
  GERGM_Object@lambda.coef <- lambda
  return(GERGM_Object)
}

