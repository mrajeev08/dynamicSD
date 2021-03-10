# Testing different arguments in simrabid --------------------------------------

# sub_cmd:=sub -t 12 -n 10 -jn test -wt 2m -md 'gdal'

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

# load in shapefile & other data
sd_shapefile <- st_read(system.file("extdata/sd_shapefile.shp", 
                                    package = "simrabid"))
load("data/sd_census_data.rda")
load("data/sd_vacc_data.rda")

# source other scriptss
source("R/sd_data.R")
source("R/utils-data.R")
source("R/summ_stats.R")

# Testing function ---------
vacc_dt <- get_sd_vacc(sd_vacc_data, sd_shapefile, origin_date = "01-Jan-2002",
                       date_fun = lubridate::dmy, units = "weeks", rollup = 4)

arg_test <- tidyr::expand_grid(sequential = c(TRUE, FALSE),
                               res_m = c(1000, 2000, 3000, NA),
                               allow_invalid = c(TRUE, FALSE),
                               leave_bounds = c(TRUE, FALSE),
                               track = c(TRUE, FALSE),
                               start_vacc = c(0, 0.2),
                               weights = c(TRUE, FALSE), # if prob get weights for each cell!
                               R0 = c(0.5, 1.1, 10),
                               iota = c(0.25, 1, NA), # if NA then pull in mock empirical incursion data
                               k = c(0.1, 1, 10),
                               break_threshold = c(0.8, 0.6),
                               vacc_dt = c("novacc", "sdvacc"), 
                               I_seeds = c(0, 3), 
                               death_rate = c(0.48, 0.65),
                               steps = c(365, 52)) # specify other arguments here (single value)


# Set up loop for testing different params

out_test <- 
  foreach(i = iter(arg_test, by = 'row'), .combine = rbind, 
          .packages = c("simrabid", "dplyr", "data.table", "sf", "raster", 
                        "magrittr")) %dopar% {
    
    if(is.na(i$res_m)) {
      admin <- TRUE 
      i$res_m <- 1000
    } else {
      admin <- FALSE
    }
    
    out <- get_sd_pops(sd_shapefile, res_m = i$res_m,
                       sd_census_data, death_rate_annual = 0.48)
    
    start_up <- setup_sim(tmax = i$steps * (5), # 5 yrs
                          rast = out$rast,
                          death_rate_annual = out$death_rate_annual,
                          birth_rate_annual = out$birth_rate_annual,
                          waning_rate_annual = 1/3,
                          block_fun = block_cells,
                          params = list(start_pop = out$start_pop),
                          step = i$steps,
                          by_admin = admin)
  
    if(i$steps == 52) {
      days_in <- 7 
      iota <- i$iota
    } else {
      days_in <- 1
      iota <- i$iota/7 # so that iotas are the same
    }
    
    param_list <- c(list(R0 = i$R0, 
                         k = i$k, 
                         iota = iota), 
                    param_defaults)
    
    
    if(i$sequential) disp <- steps_weibull else disp <- dispersal_lognorm
    
    inc_fun <- sim_incursions_pois 
    
    if(is.na(i$iota)) {
      param_list <- c(param_list, 
                      list(cell_ids = sample(start_up$cell_ids, tmax, replace = TRUE), 
                           tstep = sample(tmax, tmax, replace = TRUE)))
                      
      inc_fun <- sim_incursions_hardwired
    }
    
    if(i$weights) {
      move <- sim_movement_prob
      weights <- cell_weights(covars = list(0),
                              params = list(0),
                              start_up$ncell,
                              leave_bounds = i$leave_bounds,
                              allow_invalid = i$allow_invalid,
                              cells_block = start_up$cells_block,
                              cells_out_bounds = start_up$cells_out_bounds)
    } else {
      move <- sim_movement_continuous
      weights <- NULL
    }
    
    if(i$vacc_dt %in% "novacc") {
      vacc <- vacc_dt[0]
    } else {
      vacc <- vacc_dt
    }
    
    sim_time <- Sys.time()
    
    test <-
      tryCatch(
        expr = {
          simrabid(start_up, 
                   start_vacc = i$start_vacc, 
                   I_seeds = i$I_seeds, 
                   vacc_dt = vacc,
                   params = param_list,
                   days_in_step = days_in,
                   observe_fun = beta_detect_monthly,
                   serial_fun = serial_lognorm,
                   dispersal_fun = disp, 
                   secondary_fun = nbinom_constrained, 
                   incursion_fun = inc_fun, 
                   movement_fun = move, 
                   sequential = i$sequential, 
                   allow_invalid = i$allow_invalid,
                   leave_bounds = i$leave_bounds, 
                   max_tries = 100,
                   summary_fun = test_sim,
                   track = i$track,
                   weights = weights, 
                   row_probs = NULL,
                   coverage = FALSE,
                   break_threshold = i$break_threshold,
                   by_admin = admin) 
        },
        error = function(e){
          # write out the error & the parameter log
          err <- cbind(i, data.table(error_message = e$message))
          write.csv(err, paste0("logs/err", format(Sys.time(), "%Y%m%d"), ".csv"))
          NULL
        },
        warning = function(w){
          warn <- cbind(i, data.table(warning_message = w$message))
          write.csv(warn, paste0("logs/warn", format(Sys.time(), "%Y%m%d"), ".csv"))
          NULL
        }
      )
    
    test <- cbind(test, i, 
                  data.table(sim_time = as.numeric(Sys.time() - sim_time)))
  }

# Output results -----
write_create(out_test,
             fp("analysis/out/test/test_simrabid_opts.csv"),
             write.csv, row.names = FALSE)

# Parse these from subutil for where to put things
syncto <- "~/Documents/Projects/dynamicSD/analysis/out/"
syncfrom <- "mrajeev@della.princeton.edu:/scratch/gpfs/mrajeev/dynamicSD/analysis/out/test"

# Close out
out_session(logfile = set_up$logfile, start = set_up$start, ncores = set_up$ncores)
close_cl(cl)

print("Done:)")

