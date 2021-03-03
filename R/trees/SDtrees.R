## Transmission tree reconstruction using real Serengeti Data

##' Generate rabid data frame
##'  -----------------------------------------------------------------------------------------------
##'  1. Find cases with known dates and locations
##'  2. Exclude livestock
##'  3. Estimate uncertainty (numeric days)
##' 


##' Fit distributions
##'  -----------------------------------------------------------------------------------------------
##'  1. Serial interval for 1:n progenitors away (draw from inc + inf n times and sum and fit)
##'  2. Dispersal kernel for 1:n progenitors away (sim movement w/in SD for 1:n, and calculate 
##'     distance matrix, fit to that)
##'     

##' Guess progens and unobserved cases 1000 times
##'  -----------------------------------------------------------------------------------------------
##'  1. Generate matrices
##'  2. Get bootstrapped tree + # of unobserved cases per case
##'  

##' Look at thresholds for determining incursions (do this n x n times)
##'  ------------------------------------------------------------------------------------------------
