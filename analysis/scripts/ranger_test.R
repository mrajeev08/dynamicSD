# Ranger test

# sub_cmd:=-t 1 -n 6 -jn rtest -wt 1m -md 'gdal' -sp "analysis/scripts/ranger_test.R" -sn

# Set up on cluster ------
source("R/utils.R")
set_up <- setup_cl(mpi = FALSE)

if(!set_up$slurm) fp <- here::here else fp <- cl_safe

# packages
library(data.table)
library(abcrf)
library(sf)
library(simrabid)
library(raster)
library(dplyr)
library(magrittr)
library(foreach)
library(ranger)

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
reftab_list <-  get_reftab_list(dir = fp("analysis/out/fit"))[1:4, ]
reftl_1 <- read_reftabs(reftab_list, dir = fp("analysis/out/fit"))

# without a backend
print("this is the time without registering the cluster:\n")

print("not parallelized")
system.time(ranger::ranger(dependent.variable.name = "R0", data = reftl_1, 
                           num.threads = 1))

print("parallelized")
system.time(ranger::ranger(dependent.variable.name = "R0", data = reftl_1, 
                           num.threads = set_up$ncores))

# with the backend
print("after registering the backend:\n")
cl <- make_cl(set_up$ncores)
register_cl(cl)
print(paste("Cluster size:", cl_size(cl)))

print("not parallelized")
system.time(ranger::ranger(R0 ~ ., data = reftl_1, 
                           num.threads = 1))

print("parallelized")
system.time(ranger::ranger(R0 ~ ., data = reftl_1, 
                           num.threads = set_up$ncores))

print("also trying with the not explicit call")
system.time(ranger::ranger(dependent.variable.name = "R0", data = reftl_1, 
                           num.threads = set_up$ncores))

close_cl(cl)
print("Done!")
