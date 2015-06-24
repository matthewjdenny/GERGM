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

We are currently in the process of completing adding functionality for using node level covariates in the model and should have this functionality included in the package by the beginning of July. The model is fully functioning for networks that do not require transformation. **Note that if you are not using covariates, the network you supply to the `gergm()` function must have all edge values on the [0,1] interval.** You may do this in a number of ways, for example, we have found it works well to log and then normalize financial data so it lies on the [0,1] interval due to its heavly tailed nature.

## Example

Here is a simple working example using the `gergm( )` function: 

    library(GERGM)
    net <- matrix(runif(100),10,10)
    diag(net) <- 0
    formula <- "net ~ recip + edgeweight"  
      
    test <- gergm(formula,
              network_is_directed = TRUE,
              use_MPLE_only = FALSE,
              data_transformation = NULL,
              estimation_method = "Metropolis",
              maximum_number_of_lambda_updates = 1,
              maximum_number_of_theta_updates = 5,
              number_of_networks_to_simulate = 40000,
              thin = 1/10,
              proposal_variance = 0.5,
              exponential_weights = NULL,
              downweight_statistics_together = TRUE,
              MCMC_burnin = 10000,
              seed = 456,
              convergence_tolerance = 0.01,
              MPLE_gain_factor = 0,
              force_x_theta_update = 2,
              output_directory = getwd(),
              output_name= "Testing")

Finally you will want to check the `output_directory` which will contain a number of .pdf's which can aide in assesing model fit and in determining the statistical significance of theta parameter estimates. 

## Output

## Testing
            
So far, this package has been tested successfully on OSX 10.9.5. Please email me at <mzd5530@psu.edu> if you have success on another OS or run into any problems.
