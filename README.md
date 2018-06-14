# mfishtools

R functions for gene selection and analysis of mFISH data

mfishtools includes many functions that are used for analysis of data for the CZI SpaceTx 
project, and mostly relies on correlation-based analysis with filtering.

Install using:
```
devtools::install_github("AllenInstitute/mfishtools",auth_token="802976690281f1483c40de46d0a07e9d01a3de08")
```










#### Notes for package generation:

1) Github.  <br/>
== Get an AllenInstitute account  <br/>
== Build a blank repository  <br/>
== Add an authentication token for your personal account (check the "repo" box)  <br/>

2) Reading.  <br/>
== Rstudio package with git: https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN  <br/>
== roxygen2 https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html  <br/>
== building an R package: https://www.r-bloggers.com/building-a-package-in-rstudio-is-actually-very-easy/ (quick)  <br/>
== building an R package: http://r-pkgs.had.co.nz/ (complete)  <br/>

4) Download and install git on your computer <br/>
== Note: on Windows, you cannot install git from the desktop, it must be in a subfolder. <br/>

3) Start an R Studio project with version control (see top link in #2)  <br/>
== Link this package to the blank repo from #1  <br/>
== Note that if you want to use a network drive on windows, you need to map it  <br/>
   first (https://www.laptopmag.com/articles/map-network-drive-windows-10)  <br/>

4) Copy all of your relevant functions to the R directory  <br/>

5) Format your function annotations correctly so roxygen2 can make all the man files for you.  For example:  <br/>
```
#' Confusion matrix
#'
#' This function returns a table of the top confused clusters (assigned clusters incorrectly mapped)
#'
#' @param confusionProp confusion matrix (e.g., output from getConfusionMatrix).
#' @param count number of top confusions to show
#'
#' @return a 3 x count matrix of the top confused pairs of clusters with the three columns corresponding
#'   to mapped cluster, assigned cluster, and fraction of cells incorrectly mapped, respectively.
#'
confusion <- function()
{
  # function code here
}
```

6) Build the package  <br/>

7) Make changes and update.  Here is some useful code:  <br/>
```
roxygen2::roxygenise()  # Add the comments
formatR::tidy_dir("R",indent = getOption("formatR.indent", 2))  # Make the spacing and tabs and things consistent
```
