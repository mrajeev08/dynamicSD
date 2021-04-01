library(tibble)
library(simrabid)
library(glue)
library(readr)
par_table <- tribble(
  ~`Input (units)`, ~Distribution, ~Source,
  "Secondary case distribution (number of dogs)", 
  "Right censored negative binomial distribution with mean = $R_{0}$ and dispersion parameter $k$", 
  "$R_{0}$ and $k$ are estimated, prior centered around estimates from [@Hampson2009]", 
  "Weekly rate of introductions", "Poission distrubtion", "Estimated", 
  "Serial interval distribution (days)", 
  glue("Lognormal distribution with mean = {round(param_defaults$serial_meanlog, 2)}, and sd = {round(param_defaults$serial_sdlog, 2)}"), 
  "[@Mancyinprep]", 
  "Dispersal kernel (meters)", 
  glue("Lognormal distribution with mean = {round(param_defaults$disp_meanlog, 2)}, and sd = {round(param_defaults$disp_sdlog, 2)}"), 
  "[@Mancyinprep]", 
  "Step length distribution (meters)",
  glue("Convolved Weibull distribution with shape = {round(param_defaults$steps_shape, 2)}, and scale = {round(param_defaults$steps_scale, 2)}"), 
  "[@Mancyinprep]", 
  "Annual death rate (dogs/yr)", "1/25.8 yrs", "[@czupryna2016]", 
  "Rate of waning vaccine immunity (dogs/yr)", "1/3 yrs", "[@Hampson2009]"
)
  
write_csv(par_table, "analysis/out/par_inputs.csv")
