# Estimate parameters -----

# sub_cmd:=-t 2 -n 5 -jn test -wt 5m -md \'gdal\' -ar \'1-24\' -sn

arg <- commandArgs(trailingOnly = TRUE)

# Set up on cluster ------
source("R/utils.R")
set_up <- setup_cl(mpi = FALSE)

if(!set_up$slurm) fp <- here::here else fp <- cl_safe

# This isn't necessary for ranger but it is for the parallelization part 
# of the predict.regabcrf
cl <- make_cl(set_up$ncores)
register_cl(cl)
print(paste("Cluster size:", cl_size(cl)))

# packages
library(data.table)
library(abcrf)
library(sf)
library(simrabid)
library(raster)
library(dplyr)
library(magrittr)
library(foreach)

# scripts
source("R/utils.R")
source("R/sd_data.R")
source("R/utils-data.R")
source("R/get_observed_data.R")
source("R/summ_stats.R")
source("R/compare_mods.R")
source("R/estimate_pars.R")

# load in shapefile & other data ----
sd_shapefile <- st_read(system.file("extdata/sd_shapefile.shp", 
                                    package = "simrabid"))
load("data/sd_census_data.rda")
load("data/sd_vacc_data.rda")
incursions <- fread(fp("analysis/out/incursions.csv"))
load("data/sd_case_data.rda")

# for getting observed data (just the dummy way) ----
cand <- fread(fp("analysis/out/candidates.csv"))[1, ]
out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = cand$death_rate)
obs_data <- get_observed_data(sd_case_data, 
                              cand = cand, 
                              out = out)$obs_sstats

# set up args for array job ----
mod_ind <- as.numeric(arg[1])
cands_all <- get_name(fread(fp("analysis/out/candidates.csv")), 
                           root = TRUE)
cand_now <- unique(cands_all)[mod_ind]
reftab_list <-  data.table(get_reftab_list(dir = fp("analysis/out/abc_sims")))[root_name %in% cand_now]
reftl <- read_reftabs(reftab_list, dir = fp("analysis/out/abc_sims"))

# Param estimation ----
param_ests <- estimate_pars(reftable = reftl, 
                            par_names = c("R0", "k", "iota"), 
                            exclude = c("stopped", "sim", "break_threshold", "prop_start_pop"),
                            ncores = set_up$ncores, 
                            paral = TRUE, 
                            obs_data = obs_data, 
                            ntree = 500, 
                            predict = TRUE, 
                            predict_nsimul = 1000,
                            return_training = FALSE)

write_create(param_ests, 
             fp(paste0("analysis/out/par_ests/", cand_now, "_full.rds")),
             saveRDS)

param_ests_se <- estimate_par_se(reftable = reftl, 
                                  par_names = c("R0", "k", "iota"), 
                                  exclude = c("stopped", "sim", "break_threshold", 
                                              "prop_start_pop"),
                                  ncores = set_up$ncores, 
                                  paral = TRUE, 
                                  obs_data = obs_data, 
                                  ntree = 500, 
                                  predict = TRUE, 
                                  samp_prop = 0.75, 
                                  nsims = 3)

write_create(param_ests_se, 
             fp(paste0("analysis/out/par_ests/", cand_now, "_se.rds")),
             saveRDS)

# Parse these from subutil for where to put things
syncto <- "~/Documents/Projects/dynamicSD/analysis/out/"
syncfrom <- "mrajeev@della.princeton.edu:/scratch/gpfs/mrajeev/dynamicSD/analysis/out/par_ests"

# Close out
out_session(logfile = set_up$logfile, 
            start = set_up$start, 
            ncores = set_up$ncores)

print("Done:)")

close_cl(cl)
