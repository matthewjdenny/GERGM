Optimize_Proposal_Variance <- function(GERGM_Object,
                                       seed2,
                                       possible.stats,
                                       verbose,
                                       max_updates = 20,
                                       fine_grained_optimization = FALSE,
                                       iteration_fraction = 1/10,
                                       acceptable_bounds = 0.05,
                                       coarse_search = TRUE){

  # create a temporary GERGM object to use in proposal variance
  # optimization
  Opt_Prop_Var <- GERGM_Object
  Opt_Prop_Var@number_of_simulations <- round(max(GERGM_Object@number_of_simulations *
                                              iteration_fraction
                                            ,1000))
  Opt_Prop_Var@burnin <-  round(max(GERGM_Object@burnin * iteration_fraction
                              ,1000))
  if ((1/GERGM_Object@thin) >  (Opt_Prop_Var@number_of_simulations *
                                iteration_fraction)) {
    Opt_Prop_Var@thin <- (1/iteration_fraction)/Opt_Prop_Var@number_of_simulations
  }

  FOUND_ACCEPTABLE_PROP_VAR <- FALSE
  Acceptable_Proposal_Variance <- GERGM_Object@proposal_variance
  dampening_counter <- 1
  need_to_decrease_prop_var <- TRUE
  first_pass <- TRUE
  # find the optimal proposal variance
  while (!FOUND_ACCEPTABLE_PROP_VAR) {
    cat("--------- START HYPERPARAMETER OPTIMIZATION ---------",
        "\nSimulating",Opt_Prop_Var@number_of_simulations,
        "networks with proposal variance:", Opt_Prop_Var@proposal_variance,"\n")
    Opt_Prop_Var <- Simulate_GERGM(
      Opt_Prop_Var,
      seed1 = seed2,
      possible.stats = possible.stats,
      verbose = verbose,
      parallel = GERGM_Object@parallel_statistic_calculation)
    ar <- Opt_Prop_Var@MCMC_output$Acceptance.rate
    cat("Current acceptance rate:", ar,"\n")
    lb <- GERGM_Object@target_accept_rate - acceptable_bounds
    ub <- GERGM_Object@target_accept_rate + acceptable_bounds

    # start by doing coarse search where we change the proposal variance by
    # orders of magnitude until we over/undershoot.
    if (coarse_search) {
      # if the accept rate is too low
      if (lb > ar) {
        # if we are on our first pass, just divide the proposal variance by 10
        if (first_pass) {
          Opt_Prop_Var@proposal_variance <- Opt_Prop_Var@proposal_variance/10
          first_pass <- FALSE
        } else {
          # otherwise, if we had been increasing the proposal variance because
          # our accept rate was too high, then we know we have overshot and can
          # go to finer grained grid search.
          if (!need_to_decrease_prop_var) {
            coarse_search <- FALSE
          } else {
            # if we had not previously been increasing the porposal variance,
            # continue to make it smaller
            Opt_Prop_Var@proposal_variance <- Opt_Prop_Var@proposal_variance / 10
          }
        }
        # now deal with the case where the accept rate is too high
      } else if (ub < ar) {
        # if we are on the first pass, simply make the proposal variance 10
        # times larger
        if (first_pass) {
          Opt_Prop_Var@proposal_variance <- 10 * Opt_Prop_Var@proposal_variance
          # also make sure we note that we started out too small and need to
          # get bigger
          need_to_decrease_prop_var <- FALSE
          first_pass <- FALSE
        } else {
          # if we are not on the first pass were making it smaller, then we know
          # we have overshot and need to make it bigger.
          if (need_to_decrease_prop_var) {
            coarse_search <- FALSE
          } else {
            # otherwise keep making it bigger.
            Opt_Prop_Var@proposal_variance <- 10 * Opt_Prop_Var@proposal_variance
          }
        }
        # if we are in the sweet spot, then we do not need to keep doing coarse
        # grid search
      } else {
        coarse_search <- FALSE
      }

      # check to make sure that the proposal varince is not greater than 0.5. If
      # it is, then set it to last value and bump out of coarse search
      if (Opt_Prop_Var@proposal_variance > 0.45) {
        coarse_search <- FALSE
        Opt_Prop_Var@proposal_variance <- Opt_Prop_Var@proposal_variance / 10
      }
    } else {
      # then switch to more fine grained search once we are in the right
      # ballpark
      if (lb > ar) {
        change <- (1/(1 + dampening_counter)) * Opt_Prop_Var@proposal_variance
        Opt_Prop_Var@proposal_variance <- Opt_Prop_Var@proposal_variance - change
      } else if (ub < ar) {
        # deal with the case where making the proposal variance bigger will not help
        if (Opt_Prop_Var@proposal_variance > 0.45) {
          Acceptable_Proposal_Variance <- Opt_Prop_Var@proposal_variance
          FOUND_ACCEPTABLE_PROP_VAR <- TRUE
        } else {
          if (fine_grained_optimization) {
            # if we want to do fine grained optimization, then increase prop var by a
            # small increment
            change <- (1/(1 + dampening_counter)) * Opt_Prop_Var@proposal_variance
          } else {
            if (Opt_Prop_Var@proposal_variance > 0.25) {
              change <- (1/(1 + dampening_counter)) * (0.5 - Opt_Prop_Var@proposal_variance)
            } else {
              change <- (1/(1 + dampening_counter)) * (0.37 - Opt_Prop_Var@proposal_variance)
            }
          }
          Opt_Prop_Var@proposal_variance <- Opt_Prop_Var@proposal_variance + change
        }
      } else {
        Acceptable_Proposal_Variance <- Opt_Prop_Var@proposal_variance
        FOUND_ACCEPTABLE_PROP_VAR <- TRUE
      }

      dampening_counter <-  dampening_counter + 1
    }
    if (dampening_counter > max_updates) {
      Acceptable_Proposal_Variance <- Opt_Prop_Var@proposal_variance
      cat("Stopping optimization, more iterations will likely not improve results...\n")
      FOUND_ACCEPTABLE_PROP_VAR <- TRUE
    }
  }

  return(Acceptable_Proposal_Variance)
}
