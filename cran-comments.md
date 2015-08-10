## Test environments
* local OS X install, R 3.1.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 

* This is a new package.

## Downstream dependencies
There are no downstream dependencies as this is a new package

## Resubmission
This is a resubmission. In this version I have:

* removed header calls to math.h and cmath .
* removed the call to the ceil() function in the C++ code and replaced with a call in R code.
* included using std::xyz statements where applicable in C++ code.
