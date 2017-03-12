## ----eval=FALSE----------------------------------------------------------
#  install.packages("GERGM")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools")

## ----eval=FALSE----------------------------------------------------------
#  devtools::install_github("matthewjdenny/GERGM")

## ----eval=FALSE----------------------------------------------------------
#  library(GERGM)

## ----eval=TRUE, fig.width=6, fig.height=6, fig.align ='center'-----------
library(GERGM)
set.seed(12345)
data("lending_2005")
data("covariate_data_2005")
data("net_exports_2005")
plot_network(lending_2005) 

## ----eval=TRUE, fig.width=7, fig.height=5.5------------------------------
head(covariate_data_2005)

## ----eval=TRUE, echo=TRUE, fig.width=7, fig.height=3.5, results='hide', message=FALSE----
formula <- lending_2005 ~ edges + mutual(alpha = 0.8) + sender("log_GDP") + 
  receiver("log_GDP") + nodemix("G8", base = "No") + netcov(net_exports_2005) 

## ----eval=TRUE, echo=TRUE, fig.width=8.5, fig.height=3.5, results='hide', message=FALSE, fig.align ='center'----
test <- gergm(formula,
              covariate_data = covariate_data_2005,
	            number_of_networks_to_simulate = 40000,
	            thin = 1/100,
	            proposal_variance = 0.05,
	            MCMC_burnin = 10000,
	            seed = 456,
	            convergence_tolerance = 0.5)

## ----eval=FALSE----------------------------------------------------------
#  # Generate Estimate Plot
#  Estimate_Plot(test)
#  # Generate GOF Plot
#  GOF(test)
#  # Generate Trace Plot
#  Trace_Plot(test)

## ----eval=TRUE, echo=TRUE, fig.width=6.5, fig.height=3, results='hide', message=FALSE, fig.align ='center'----
Estimate_Plot(test,
              coefficients_to_plot = "both",
              coefficient_names = c("Mutual Dyads",
                                    "log(GDP) Sender",
                                    "log(GDP) Receiver",
                                    "Non-G8 Sender, G8 Receiver",
                                    "G8 Sender, Non-G8 Receiver",
                                    "G8 Sender, G8 Receiver",
                                    "intercept",
                                    "Normalized Net Exports",
                                    "Dispersion Parameter"),
              leave_out_coefficients = "intercept")

## ----eval=TRUE, echo=TRUE, fig.width=5, fig.height=5, results='hide', message=FALSE, fig.align ='center'----
# Generate Hysteresis plots for all structural parameter estimates
hysteresis_results <- hysteresis(test,
                                 networks_to_simulate = 1000,
                                 burnin = 300,
                                 range = 8,
                                 steps = 20,
                                 simulation_method = "Metropolis",
                                 proposal_variance = 0.05)

## ----eval=TRUE, echo=TRUE, results='hide', message=FALSE-----------------
test2 <- conditional_edge_prediction(
  GERGM_Object = test,
  number_of_networks_to_simulate = 100,
  thin = 1,
  proposal_variance = 0.05,
  MCMC_burnin = 100,
  seed = 123)

## ----eval=TRUE-----------------------------------------------------------
MSE_results <- conditional_edge_prediction_MSE(test2)

## ----eval=FALSE----------------------------------------------------------
#    set.seed(12345)
#    # Function to generating a random positive-definite matrix with user-specified
#    # positive eigenvalues. If eigenvalues are not specified, they are generated
#    # from a uniform distribution.
#    Posdef <- function (n, ev = runif(n, 0, 10)) {
#      Z <- matrix(ncol=n, rnorm(n^2))
#      decomp <- qr(Z)
#      Q <- qr.Q(decomp)
#      R <- qr.R(decomp)
#      d <- diag(R)
#      ph <- d / abs(d)
#      O <- Q %*% diag(ph)
#      Z <- t(O) %*% diag(ev) %*% O
#      return(Z)
#    }
#  
#    # Generate eignevalues
#    x <- rnorm(10)
#    # generate a positive definite matrix
#    pdmat <- Posdef(n = 10)
#    # transform to correlations
#    correlations <- pdmat / max(abs(pdmat))
#    diag(correlations) <- 1
#    net <- (correlations + t(correlations)) / 2
#  
#    # add in node names
#    colnames(net) <- rownames(net) <- letters[1:10]
#  
#    # correlation GERGM specification
#    formula <- net ~ edges + ttriads
#  
#    # model should run in under a minute
#    test <- gergm(formula,
#                  estimation_method = "Metropolis",
#                  number_of_networks_to_simulate = 100000,
#                  thin = 1/100,
#                  proposal_variance = 0.2,
#                  MCMC_burnin = 100000,
#                  seed = 456,
#                  convergence_tolerance = 0.5,
#                  beta_correlation_model = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  formula <- net ~  mutual(0.8) + ttriads(0.8) + out2stars(0.8) +
#    sender("log_GDP") + netcov(net_exports) +
#    receiver("log_GDP") + nodemix("G8", base = "No")
#  
#  
#  result <- gergm(formula,
#                covariate_data = covariate_data_2005,
#                number_of_networks_to_simulate = 400000,
#                thin = 1/100,
#                proposal_variance = 0.05,
#                MCMC_burnin = 200000,
#                seed = 456,
#                convergence_tolerance = 0.8,
#                hyperparameter_optimization = TRUE,
#                target_accept_rate = 0.25,
#                weighted_MPLE = TRUE,
#                theta_grid_optimization_list = list(grid_steps = 2,
#                                                    step_size = 0.1,
#                                                    cores = 30,
#                                                    iteration_fraction = 1))

