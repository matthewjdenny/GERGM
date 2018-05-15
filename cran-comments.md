## Test environments
* Local OS X install, R 3.5.0 and devel
* Winbuilder (devel and release)
* Ubuntu precise 14.04.5 LTS (travis CI), R 3.5.0

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs on my installation of Windows, OS X, or Linux. This package also fixes the issues found in the CRAN checks. To the best of my knowledge, the NOTE about the size of the libs directory in the installed package is due to its reliance on the BH package, which is important for the estimation routines it implements, and does not appear to be unique to this package.

## Downstream dependencies
The package update does not cause any issues for the one package that depends on it, to the best of my knowledge. 


