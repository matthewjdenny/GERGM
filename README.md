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

    install.packages( pkgs = c("BH","RcppArmadillo"), dependencies = TRUE)

If all went well, check out the `?GERGM` help file to see a full working example with info on how the data should look. 

## Basic Useage

## Example


## Output

## Testing
            
So far, this package has been tested successfully on OSX 10.9.5. Please email me at <mzd5530@psu.edu> if you have success on another OS or run into any problems.
