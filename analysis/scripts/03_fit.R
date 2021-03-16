# Candidate models (to benchmark) ---------

# sub_cmd:=-t 1 -n 21 -jn fit -wt 1m -md 'gdal' -ar '1-4' -cmd '100'

arg <- commandArgs(trailingOnly = TRUE)

# Set up on cluster ------
source("R/utils.R")
set_up <- setup_cl(mpi = TRUE)

cl <- make_cl(set_up$ncores)
register_cl(cl)
print(paste("Cluster size:", cl_size(cl)))

if(!set_up$slurm) fp <- here::here else fp <- cl_safe

# Dependencies
library(raster)
library(data.table)
library(sf)
library(tidyr)
library(dplyr)
library(magrittr)
library(simrabid) # remotes::install_github("mrajeev08/simrabid")
library(foreach)
library(iterators)
library(doRNG)
library(lubridate)

# set up args
mod_ind <- as.numeric(arg[1])
cand <- fread(fp("analysis/out/fit/candidates.csv"))[mod_ind, ]
nsims <- as.numeric(arg[2])

# load in shapefile & other data
sd_shapefile <- st_read(system.file("extdata/sd_shapefile.shp", 
                                    package = "simrabid"))
load("data/sd_census_data.rda")
load("data/sd_vacc_data.rda")
load("data/incursions.rda")
load("data/sd_case_data.rda")

# source other scriptss
source("R/sd_data.R")
source("R/get_observed_data.R")
source("R/utils-data.R")
source("R/summ_stats.R")

# Testing function ---------
vacc_dt <- get_sd_vacc(sd_vacc_data, sd_shapefile, origin_date = "01-Jan-2002",
                       date_fun = lubridate::dmy, units = "weeks", rollup = 4)

# Get the startup (wrap everything below this in a function)
# pass cand & summary stats & list of param values (this is where it's tricky)
# maybe check against default params? if not there then pull from 
# list that gets passed (should either be length 1 or length nsims)
out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = cand$death_rate)

start_up <- setup_sim(start_date = "2002-01-01", 
                      apprx_end_date = "2020-12-31", 
                      days_in_step = 7, 
                      rast = out$rast,
                      death_rate_annual = out$death_rate_annual,
                      birth_rate_annual = out$birth_rate_annual,
                      waning_rate_annual = 1/3,
                      block_fun = block_cells,
                      params = list(start_pop = out$start_pop),
                      by_admin = cand$by_admin)

if(!cand$estincs) {
  # add cell ids & timesteps of incursion to parameter defaults
  incursions %<>%
    mutate(cell_id = raster::cellFromXY(out$rast,
                                        cbind(x_coord, y_coord)), 
           tstep = get_timestep(date, 
                                origin_date = "2002-01-01",
                                date_fun = lubridate::ymd, 
                                units = "weeks"))
  mean_per_week <- mean(tabulate(incursions$tstep))
  sim_inc_steps <- seq(ceiling(max(incursions$tstep)), start_up$tmax)
  sim_inc_n <- rpois(length(sim_inc_steps), mean_per_week)
  tstep_add <- rep(sim_inc_steps, sim_inc_n)
  cell_ids_add <- sample(start_up$cell_ids, length(tstep_add))
  
  param_defaults <- c(param_defaults, 
                      list(cell_ids = c(incursions$cell_id, cell_ids_add), 
                           tstep = c(incursions$tstep, tstep_add)))
}

if(cand$weights) {
  weights <- cell_weights(covars = list(0),
                          params = list(0),
                          start_up$ncell,
                          leave_bounds = TRUE,
                          allow_invalid = FALSE,
                          cells_block = start_up$cells_block,
                          cells_out_bounds = start_up$cells_out_bounds)
} else {
  weights <- NULL
}

# params (do this parallel, spit out i + systemtime + wrap in trycatch!
disp_fn <- ifelse(cand$sequential, simrabid::steps_weibull, simrabid::dispersal_lognorm)
inc_fn <- ifelse(cand$estincs, simrabid::sim_incursions_pois, simrabid::sim_incursions_hardwired)
move_fn <- ifelse(cand$weights, simrabid::sim_movement_prob, simrabid::sim_movement_continuous)

# Set up priors
priors <- list(R0 = function(n) exp(rnorm(n, mean = 0.1, sd = 0.2)), # centered around 1.1
               iota = function(n) exp(rnorm(n, mean = 0.5, sd = 0.5)), # centered around 1.5
               k = function(n) exp(rnorm(n, mean = 0.5, sd = 0.5))) # uniform

if(!cand$estincs) priors$k <- function(n) rep(0, n)

# get observed data
obs_data <- get_observed_data(sd_case_data, 
                              mod_specs = cand, 
                              mod_start = out)$obs_data

out_sims <-
  foreach(j = seq_len(nsims), .combine = rbind, 
        .packages = c("simrabid", "dplyr", 
                      "data.table", "sf", "raster", 
                      "magrittr"), 
        .options.RNG = cand$seed) %dorng% {
        
       # draw a val 
       R0_val <- priors$R0(1)
       k_val <- priors$k(1)
       iota_val <- priors$iota(1) 
          
        out <- 
          tryCatch(
            expr = {
              simstats <- simrabid(start_up, 
                                   start_vacc = cand$start_vacc, 
                                   I_seeds = cand$I_seeds, 
                                   vacc_dt = vacc_dt,
                                   params = c(R0 = R0_val, k = k_val, 
                                              iota = iota_val, 
                                              param_defaults),
                                   days_in_step = cand$days_in_step,
                                   observe_fun = beta_detect_monthly,
                                   serial_fun = serial_lognorm,
                                   dispersal_fun = disp_fn, 
                                   secondary_fun = nbinom_constrained, 
                                   incursion_fun = inc_fn, 
                                   movement_fun = move_fn,
                                   sequential = cand$sequential, 
                                   allow_invalid = cand$allow_invalid,
                                   leave_bounds = cand$leave_bounds, 
                                   max_tries = 100,
                                   summary_fun = inc_stats, # use inc_stats after testing!
                                   track = cand$track,
                                   weights = weights, 
                                   row_probs = NULL,
                                   coverage = FALSE,
                                   break_threshold = cand$break_threshold,
                                   by_admin = cand$by_admin,
                                   extra_pars = list(obs_data = obs_data))
              
              out <- data.table(simstats, R0 = R0_val, 
                                k = k_val, 
                                iota = iota_val)
            },
            error = function(e){
              # write out the error & the parameter log
              err <- cbind(R0 = R0_val, k = k_val, 
                           iota = iota_val, cand,
                           e$message)
              write.csv(err, 
                        paste0("logs/err", 
                               format(Sys.time(), "%Y%m%d_%H%M%S"), 
                               ".csv"),
                        row.names = FALSE)
              NULL
            }
          )
  } 

file_out <- paste0("analysis/out/fit/", get_name(cand), ".csv")


write_create(out_sims,
             fp(file_out),
             data.table::fwrite)


# Parse these from subutil for where to put things
syncto <- "~/Documents/Projects/dynamicSD/analysis/out/"
syncfrom <- "mrajeev@della.princeton.edu:/scratch/gpfs/mrajeev/dynamicSD/analysis/out/fit"

# Close out
out_session(logfile = set_up$logfile, start = set_up$start, ncores = set_up$ncores)
close_cl(cl)

print("Done:)")