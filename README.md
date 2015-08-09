# GERGM
An R package to estimate Generalized Exponential Random Graph Models

NOTE: **This package is still under development and is very much a work in progress. Please do not use in any publication without express consent of the authors. PLEASE REPORT ANY BUGS OR ERRORS TO <mzd5530@psu.edu>**. 

## Model Overview 

An R package which implements the Generalized Exponential Random Graph Model (GERGM) with an extension to estimation via Metropolis Hastings. The relevant papers detailing the model can be found at the links below:

* Bruce A. Desmarais, and Skyler J. Cranmer,  (2012). "Statistical inference for valued-edge networks: the generalized exponential random graph model". PloS One. [[Available Here](http://dx.plos.org/10.1371/journal.pone.0030136)]
* James D. Wilson, Matthew J. Denny, Shankar Bhamidi, Skyler Cranmer, and Bruce Desmarais (2015). "Stochastic Weighted Graphs: Flexible Model Specification and Simulation". [[Available Here](http://arxiv.org/abs/1505.04015)]

## Installation

### Requirements for using C++ code with R

See the **Requirements for using C++ code with R** section in the following tutorial: [Using C++ and R code Together with Rcpp](http://www.mjdenny.com/Rcpp_Intro.html).

### Installing The Package
  
To install this package from Github, you will need to Hadley Wickham's devtools package installed.

    install.packages("devtools")
    library("devtools")
    
Now we can install from Github using the following line:

    devtools::install_github("matthewjdenny/GERGM")

I have  had success installing with R 3.2.0+ installed but if you do not have the latest version of R installed, or run into some install errors (please email if you do), it should work as long as you install the dependencies first with the following block of code:

    install.packages( pkgs = c("BH","RcppArmadillo","ggplot2","methods"), dependencies = TRUE)

Once the `GERGM` package is installed, you may access its functionality as you would any other package by calling:

    library(GERGM)

If all went well, check out the `?GERGM` help file to see a full working example with info on how the data should look. 

## Basic Useage

To use this package, first load in the network you wish to use as a (square) matrix, following the example provided below. You may then use the gergm() function to estimate a model using any combination of the following statistics: `out2star(alpha = 1)`, `in2star(alpha = 1)`, `ctriads(alpha = 1)`, `recip(alpha = 1)`, `ttriads(alpha = 1)`, `edges(alpha = 1)`, `absdiff(covariate = "MyCov")`, `edgecov(covariate = "MyCov")`, `sender(covariate = "MyCov")`, `reciever(covariate = "MyCov")`, `nodefactor(covariate, base = "MyBase")`, `netcov(network_covariate)`. To use exponential downweighting for any of the network level terms, simply specify a value for alpha less than 1. The `gergm()` function will provide all of the estimation and diagnostic functionality and the parameters of this function can be querried by typing `?gergm` into the R console. You may also generate diagnostic plots using a GERGM Object returned by the `gergm()` function by using any of the following functions: `Estimate_Plot()`, `GOF()`, `Trace_Plot()`.

## Examples

Here are two simple working examples using the `gergm()` function: 

    library(GERGM)
    ########################### 1. No Covariates #############################
    # Preparing an unbounded network without covariates for gergm estimation #
    net <- matrix(rnorm(100,0,20),10,10)
    colnames(net) <- rownames(net) <- letters[1:10]
    formula <- net ~ recip + edges  
      
    test <- gergm(formula,
    	          normalization_type = "division",
    	          network_is_directed = TRUE,
    	          use_MPLE_only = FALSE,
    	          estimation_method = "Metropolis",
    	          maximum_number_of_lambda_updates = 1,
    	          maximum_number_of_theta_updates = 5,
    	          number_of_networks_to_simulate = 40000,
    	          thin = 1/10,
    	          proposal_variance = 0.5,
    	          downweight_statistics_together = TRUE,
    	          MCMC_burnin = 10000,
    	          seed = 456,
    	          convergence_tolerance = 0.01,
    	          MPLE_gain_factor = 0,
    	          force_x_theta_update = 4)
      
    ########################### 2. Covariates #############################
    # Preparing an unbounded network with covariates for gergm estimation #
    net <- matrix(runif(100,0,1),10,10)
    colnames(net) <- rownames(net) <- letters[1:10]
    node_level_covariates <- data.frame(Age = c(25,30,34,27,36,39,27,28,35,40),
    	                                Height = c(70,70,67,58,65,67,64,74,76,80),
    	                                Type = c("A","B","B","A","A","A","B","B","C","C"))
    rownames(node_level_covariates) <- letters[1:10]
    network_covariate <- net + matrix(rnorm(100,0,.5),10,10)
    formula <- net ~ recip + edges + sender("Age") + 
    netcov("network_covariate") + nodefactor("Type",base = "A")  
       
    test <- gergm(formula,
    	          covariate_data = node_level_covariates,
    	          network_is_directed = TRUE,
    	          use_MPLE_only = FALSE,
    	          estimation_method = "Metropolis",
    	          maximum_number_of_lambda_updates = 5,
    	          maximum_number_of_theta_updates = 5,
    	          number_of_networks_to_simulate = 100000,
    	          thin = 1/10,
    	          proposal_variance = 0.5,
    	          downweight_statistics_together = TRUE,
    	          MCMC_burnin = 50000,
    	          seed = 456,
    	          convergence_tolerance = 0.01,
    	          MPLE_gain_factor = 0,
    	          force_x_theta_update = 2)
      
    # Generate Estimate Plot
    Estimate_Plot(test)
    # Generate GOF Plot
    GOF(test)
    # Generate Trace Plot
    Trace_Plot(test)

Finally, if you specified an `output_directory` and `output_name`, you will want to check the `output_directory` which will contain a number of .pdf's which can aide in assesing model fit and in determining the statistical significance of theta parameter estimates. 

### Output

If `output_name` is specified in the `gergm()` function, then five files will be automatically generated and saved to the `output_directory`. The example file names provided below are for `output_name = "Test"`:

* **"Test_GOF.pdf"**  -- Normalized, standardized Goodness Of Fit plot comparing statistics of networks simulated from final parameter estimates to the observed network.
* **"Test_Parameter_Estimates.pdf"** -- Theta parameter estimates with 90 and 90% confidence intervals.
* **"Test_GERGM_Object.Rdata"** -- The GERGM object returned by the `gergm()` functionafter estimation.
* **"Test_Estimation_Log.txt"** -- A log of all console output generated by the `gergm()` function.
* **"Test_Trace_Plot.pdf"** -- Trace plot from last round of network simulations used to generate GOF plot, useful for diagnosing an acceotance rate that is too low.


## Testing
            
So far, this package has been tested successfully on OSX 10.9.5 and Windows 7. Please email me at <mzd5530@psu.edu> if you have success on another OS or run into any problems.
