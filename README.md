# EPAbvs

This folder contains the code for the Rcpp package EPAbvs for the manuscript:
Functional Concurrent Regression Mixture Models Using Spiked Ewens-Pitman Attraction Priors
by Mingrui Liang, Matthew D. Koslovsky and Marina Vannucci

The main functions operate in C++ via the R package Rcpp. 
These functions can be sourced by a set of wrapper functions that enable easy implementation of the code in the R environment. 
Various functions are available that produce, summarize, and plot the results for inference.  
This package relies on various R packages that need to be installed in the R environment before running. 
To install, use the install.packages("") command for the following packages:  

  Rcpp   
  RcppArmadillo  
  mcclust  
  spikeSlabGAM  
  MCMCpack  
  mvtnorm  
  ggplot2
and also 
  devtools

To install the package, download the ‘EPAbvs’ folder and set the working directory in Rstudio to its path, then run  
  library( devtools )  
  devtools::build(vignettes = F)  
  devtools::install()  
  
A worked example and technical details are provided in vignettes/VignetteEPABVS.html  
