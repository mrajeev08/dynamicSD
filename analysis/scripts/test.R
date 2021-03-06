library(raster)
library(data.table)
library(sf)
library(tidyr)
library(dplyr)
library(magrittr)
devtools::load_all() # with this data gets lazy loaded!
# library(simrabid)


# R0 prior (exp of uniform -0.5 to 0.5 translates to 0.6 - 1.6 or norm centered around 1)
range(exp(rnorm(10000, mean = 0, sd = 0.2)))
# k prior (centered around 1?)
range(exp(rnorm(1000, mean = 0, sd = 2)))
# incursions
range(exp(rnorm(1000, mean = 0, sd = 1)))


# candidates
test_df <- expand.grid(leave_bounds = c(TRUE, FALSE),
                       allow_invalid = c(TRUE, FALSE),
                       movement = c("prob", "continuous"),
                       res_m = c(1000, 2000, 4000, "vill"))

# 1000 sims for each across the priors

# output = timings + distance kernel + prop of movements that were to invalid/out-of-bounds

# for a specific estimate (i.e. what I think it is) do 1000 sims and get timeseries

# Testing function
res_test <- c(1000, 2000, 4000, "admin")
arg_test <- expand_grid(sequential = c(TRUE, FALSE),
                        allow_invalid = c(TRUE, FALSE),
                        leave_bounds = c(TRUE, FALSE),
                        track = c(TRUE, FALSE),
                        move = c(sim_movement_continuous, sim_movement_prob), # if prob get weights for each cell!
                        R0 = c(0.5, 1.1, 1.5),
                        iota = c(0.25, 1, 5, NA), # if NA then pull in mock empirical incursion data
                        k = c(0.1, 1, 10),
                        break_threshold = c(0.25, 0.5, 0.8)) # specify other arguments here (single value)

# summary function to compare time series & spatial distribution of cases

# set up the loop here (starting with the resolution decision)
# 1. get pop aggregated to the correct resolution
out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = 0.48)

# vacc data remains constant throughout this
vacc_dt <- get_sd_vacc(sd_vacc_data, sd_shapefile, origin_date = "01-Jan-2002",
                       date_fun = lubridate::dmy, units = "weeks", rollup = 4)

# 2. set-up simulation framework
start_up <- setup_sim(tmax = 52* (2020 - 2002),
                      rast = out$rast,
                      death_rate_annual = out$death_rate_annual,
                      birth_rate_annual = out$birth_rate_annual,
                      waning_rate_annual = 1/3,
                      block_fun = block_cells,
                      params = list(start_pop = out$start_pop),
                      step = 52,
                      by_admin = FALSE)
weights <- cell_weights(covars = list(0),
                        params = list(0),
                        start_up$ncell,
                        leave_bounds = TRUE,
                        allow_invalid = FALSE,
                        cells_block = start_up$cells_block,
                        cells_out_bounds = start_up$cells_out_bounds)
test <-
  tryCatch(
    expr = {
      simrabid(start_up, start_vacc = 0.2, I_seeds = 3, vacc_dt,
               params = c(list(R0 = 1.2, k = 1, iota = 1),
                          param_defaults),
               days_in_step = 7,
               observe_fun = beta_detect_monthly,
               serial_fun = serial_lognorm,
               dispersal_fun = steps_weibull,
               secondary_fun = nbinom_constrained,
               incursion_fun = sim_incursions_pois,
               movement_fun = sim_movement_prob,
               sequential = TRUE, allow_invalid = TRUE,
               leave_bounds = TRUE, max_tries = 100,
               summary_fun = check_sim,
               track = TRUE,
               weights = weights,
               row_probs = NULL,
               coverage = FALSE,
               break_threshold = 0.8,
               by_admin = FALSE)
    },
    error = function(e){
      message('Caught an error!') # write out the error & the parameter log
      print(e)
      NULL
    },
    warning = function(w){
      message('Caught an warning!')  # write out the warning & the parameter log
      print(w)
      NULL
    }
  )

