## Test environments
* local OS X install, R 3.2.1
* win-builder (devel and release)
* travis CI, R 3.2.1

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 

* This is a new package.

## Downstream dependencies
There are no downstream dependencies as this is a new package

## Resubmission
This is a resubmission. In this version I have:

* changed from using "boost::uniform_real_distribution" to using "boost::uniform_01" which removes assert-fail WARNING when compiling under linux (ubuntu) using gcc. This shoudl address the issue noted by Dr. Hornik. 
