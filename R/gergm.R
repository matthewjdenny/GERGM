
# main function
gergm <- function(formula_object,
                  directed = c(TRUE, FALSE),
                  MPLE.only = c(FALSE, TRUE),
                  transform.data = NULL,
                  method = c("Gibbs", "Metropolis"),
                  max.num.iterations = 10,
                  mc.num.iterations = 100,
                  nsim = 500,
                  thin = 1,
                  shape.parameter = 1,
                  exponential_weights = NULL,
                  together = 0,
                  MCMC.burnin = 100,
                  seed = 123,
                  tolerance = 0.01,
                  gain.factor = 0){


  #1. Create GERGM object from network


  #2. Estimate GERGM

  Estimate_GERGM

  #3. Perform degeneracy diagnostics and create GOF plots


  #4. Return GERGM object


}
