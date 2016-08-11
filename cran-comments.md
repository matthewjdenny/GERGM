## Test environments
* local OS X install, R 3.3.1
* win-builder (devel and release)
* travis CI, R 3.3.1

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs on Windows, OS X, or Linux. However, on CRAN Checks, r-devel-linux-x86_64-fedora-clang produces a warning. I have investigated this warning and it derives from a few lines in the source code of the RcppParalell package (which my package depends on), to the best of my knowledge.

## Downstream dependencies
The package update does not cause any issues for the one package that depends on it, to the best of my knowledge. 


