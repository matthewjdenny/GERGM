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

## ----eval=TRUE, echo=TRUE, results='hide', message=FALSE----
test2 <- conditional_edge_prediction(
  GERGM_Object = test,
  number_of_networks_to_simulate = 100,
  thin = 1,
  proposal_variance = 0.05,
  MCMC_burnin = 100,
  seed = 123)

## ----eval=TRUE---------------------------------------------
MSE_results <- conditional_edge_prediction_MSE(test2)

## ----eval=FALSE--------------------------------------------
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

