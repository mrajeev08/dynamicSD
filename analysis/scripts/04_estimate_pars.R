# Estimate parameters

estimate_pars <- function(reftable, 
                          par_names = c("R0", "k", "iota"), 
                          exclude = c("stopped", "sim", "break_threshold"),
                          ncores, 
                          paral, 
                          obs_data = NULL, 
                          sampsize = nrow(reftable), 
                          ntree = 500, 
                          predict = TRUE, 
                          return_training = FALSE) {
  
  # First get reftable set up right
  
  # Then loop through and for each parameter
  # do the estimation
  # get out the var importance, the oob errs, the weights 
  # pass the par name to the df
  
  
  
}
