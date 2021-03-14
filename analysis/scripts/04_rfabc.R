# Random Forest ABC: Parameter estimation

# packages
library(data.table)
library(abcrf)
library(sf)
library(simrabid)
library(raster)
library(dplyr)
library(magrittr)

# scripts
source("R/utils.R")
source("R/sd_data.R")
source("R/utils-data.R")
source("R/get_observed_data.R")
source("R/summ_stats.R")
source("R/density_mod.R")

# load in shapefile & other data
sd_shapefile <- st_read(system.file("extdata/sd_shapefile.shp", 
                                    package = "simrabid"))
load("data/sd_census_data.rda")
load("data/sd_vacc_data.rda")
load("data/incursions.rda")
load("data/sd_case_data.rda")

# model of interest
cand <- fread("analysis/out/fit/candidates.csv")[11, ]
reftab_name <- get_name(cand)
reftable_all <- fread(paste0("analysis/out/fit/", reftab_name, ".csv"))

# for getting observed data
out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = cand$death_rate)

# get observed data
obs_data <- get_observed_data(sd_case_data, 
                              mod_specs = cand, 
                              mod_start = out)$obs_sstats

# par estimation
reftable <- reftable_all[sample(nrow(reftable_all), 1e4)]

model_rf_iota <- regAbcrf(iota ~ ., reftable[, !c("R0", "k")], 
                          paral = TRUE, ncores = 3)
plot(model_rf_iota)
predict(model_rf_iota, obs_data, reftable[, !c("R0", "k")])
out <- weights(model_rf_iota, obs_data, reftable[, !c("R0", "k")])

# Try on cluster w & without detect_cores & doMPI backend (check how the paral works in rfabv)
model_rf_R0 <- regAbcrf(R0 ~ ., reftable[, !c("iota", "k")], 
                        paral = TRUE, ncores = 3) # this works
plot(model_rf_R0)
predict(model_rf_R0, obs_data, reftable[, !c("iota", "k")], paral = TRUE, 
        ncores = 3)
out <- weights(model_rf_R0, obs_data, reftable[, !c("iota", "k")], paral = TRUE, 
               ncores = 3)

model_rf_k <- regAbcrf(k ~ ., reftable[, !c("R0", "iota")], paral = TRUE, 
                       ncores = 3)
plot(model_rf_k)
predict(model_rf_k, obs_data, reftable[, !c("R0", "iota")], paral = TRUE)
out2 <- weights(model_rf_k, obs_data, reftable[, !c("R0", "iota")], ncores = 3)

covRegAbcrf(model_rf_iota, model_rf_R0, obs_data, reftable[, !c("R0", "k")],
            reftable[, !c("iota", "k")], paral = TRUE, ncores = 3)

covRegAbcrf(model_rf_k, model_rf_R0, obs_data, reftable[, !c("R0", "iota")],
            reftable[, !c("iota", "k")], paral = TRUE, ncores = 3)

covRegAbcrf(model_rf_k, model_rf_iota, obs_data, reftable[, !c("R0", "iota")],
            reftable[, !c("R0", "k")], paral = TRUE, ncores = 3)


# model comparison
# model of interest
cand <- fread("analysis/out/fit/candidates.csv")[9, ]
reftab_name <- get_name(cand)
reftable_all <- fread(paste0("analysis/out/fit/", reftab_name, ".csv"))
reftable <- reftable_all[sample(nrow(reftable_all), 1e4)]
reftable <- reftable[, !c("R0", "k", "iota")][, modindex := factor(1)]

cand2 <- fread("analysis/out/fit/candidates.csv")[11, ]
reftab_name2 <- get_name(cand2)
reftable_all2 <- fread(paste0("analysis/out/fit/", reftab_name2, ".csv"))
reftable2 <- reftable_all2[sample(nrow(reftable_all2), 1e4)]
reftable2 <- reftable2[, !c("R0", "k", "iota")][, modindex := factor(2)]

cand3 <- fread("analysis/out/fit/candidates.csv")[13, ]
reftab_name3 <- get_name(cand3)
reftable_all3 <- fread(paste0("analysis/out/fit/", reftab_name3, ".csv"))
reftable3 <- reftable_all3[sample(nrow(reftable_all3), 1e4)]
reftable3 <- reftable3[, !c("R0", "k", "iota")][, modindex := factor(3)]

cand4 <- fread("analysis/out/fit/candidates.csv")[15, ]
reftab_name4 <- get_name(cand4)
reftable_all4 <- fread(paste0("analysis/out/fit/", reftab_name4, ".csv"))
reftable4 <- reftable_all4[sample(nrow(reftable_all4), 1e4)]
reftable4 <- reftable4[, !c("R0", "k", "iota")][, modindex := factor(4)]

reftab <- rbind(reftable, reftable2, reftable3, reftable4)

mod_comp <- abcrf(modindex~., data = reftab, 
                  group = list("1", "2", "3", "4"), 
                  paral = TRUE, ncores = 3)
predict(mod_comp, obs_data, reftab, paral = TRUE, ncores = 3)
plot(mod_comp, reftab, paral = TRUE, ncores = 3)
plot(mod_comp, reftab, obs = obs_data, paral = TRUE, ncores = 3)
