# Random Forest ABC: Model choice -------

# sub_cmd:=-t 6 -n 12 -jn comp -wt 5m -md \"gdal\" -sn -@

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

# load in shapefile & other data
sd_shapefile <- st_read(system.file("extdata/sd_shapefile.shp", 
                                    package = "simrabid"))
load("data/sd_census_data.rda")
load("data/sd_vacc_data.rda")
load("data/incursions.rda")
load("data/sd_case_data.rda")

# for getting observed data
cand <- fread(fp("analysis/out/fit/candidates.csv"))[1, ]
out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = cand$death_rate)
obs_data <- get_observed_data(sd_case_data, 
                              cand = cand, 
                              out = out)$obs_sstats

# Model comparison -----------
reftab_list <-  get_reftab_list(dir = fp("analysis/out/fit"))
reftl <- read_reftabs(reftab_list, dir = fp("analysis/out/fit"))

mod_comp <- compare_mods(reftable = reftl, 
                         par_names = c("R0", "k", "iota"), 
                         exclude = c("stopped", "sim", "break_threshold"),
                         obs_data = obs_data, 
                         ntree = 500, 
                         ncores = set_up$ncores, 
                         paral = TRUE,
                         predict = TRUE, 
                         return_training = FALSE)

saveRDS(mod_comp, fp("analysis/out/mod_comp_full.rds"))

mod_comp_se <- compare_mod_se(reftable = reftl, 
                              par_names = c("R0", "k", "iota"), 
                              exclude = c("stopped", "sim", "break_threshold"),
                              obs_data = obs_data, 
                              samp_prop = 0.75, 
                              nsims = 3, 
                              ntree = 500, 
                              ncores = set_up$ncores, 
                              paral = TRUE, 
                              predict = TRUE,
                              return_training = FALSE) 

saveRDS(mod_comp_se, fp("analysis/out/mod_comp_se.rds"))

# Parse these from subutil for where to put things
syncto <- "~/Documents/Projects/dynamicSD/analysis/"
syncfrom <- "mrajeev@della.princeton.edu:/scratch/gpfs/mrajeev/dynamicSD/analysis/out"

# Close out
out_session(logfile = set_up$logfile, 
            start = set_up$start, 
            ncores = set_up$ncores)

print("Done:)")

close_cl(cl)

