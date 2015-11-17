## Resubmission
This is a resubmission. In this version I have:

* Prevented a long-running example from running in R CMD check as requested by Dr. Hornik. 

* Attempted to address the non-portabiliy of the PKG_LIBS call to -llapack by replacing it.

* In response to Dr. Ripley's comment about the example code failing without long doubles, I
was able to replicate this issue and beleive it is due to the particular specification of 
the model in the example. While this code is no longer being run in R CMD check, this error
is not uninformative, as it indicates to the user that the model they are estimating has 
become degenerate. 

## Test environments
* local OS X install, R 3.2.2
* win-builder (devel and release)
* travis CI, R 3.2.2

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 

* I get a note when submitting using devtools because I have updated my email address.

## Downstream dependencies
There are no downstream dependencies as this is a new package


