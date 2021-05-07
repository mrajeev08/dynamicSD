# Do transmission tree reconstruction & generate priors -----

source("R/utils.R")

# Packages
library(here)
library(treerabid) # devtools::install_github("mrajeev08/treerabid")
library(data.table)
library(lubridate)
library(dplyr)
library(lubridate)
library(magrittr)
library(foreach)
library(iterators)
library(doRNG)
library(simrabid)

# Data ----
load(file = "data/sd_case_data_trees.rda")

# clean up (no cases with NA location or time & filter to start/end dates)
sd_case_data_trees %>%
  mutate(symptoms_started = dmy(symptoms_started)) %>%
  dplyr::filter(!is.na(symptoms_started), 
                !is.na(utm_easting), 
                !is.na(utm_easting), 
                symptoms_started >= ymd("2002-01-01"),
                symptoms_started <= ymd("2020-12-31")) %>%
  # get uncertainty in days
  mutate(days_uncertain = case_when(symptoms_started_accuracy == "+/- 14 days" ~ 14L, 
                                    symptoms_started_accuracy == "+/- 7 days" ~ 7L,
                                    symptoms_started_accuracy == "+/- 28 days" ~ 28L, 
                                    symptoms_started_accuracy == "0" ~ 0L, 
                                    TRUE ~ 0L)) -> case_dt 

# Get reconstructed tree & incursions ----
cl <- parallel::makeCluster(parallel::detectCores() - 1)
doParallel::registerDoParallel(cl)
ttrees <- 
  boot_trees(id_case = case_dt$id,
             id_biter = case_dt$biter_id, 
             x_coord = case_dt$utm_easting,
             y_coord = case_dt$utm_northing,
             owned = FALSE, 
             date_symptoms = case_dt$symptoms_started, # needs to be in a date class
             days_uncertain = case_dt$days_uncertain,
             use_known_source = TRUE,
             prune = TRUE,
             si_fun = treerabid::si_lnorm1,
             dist_fun = treerabid::dist_lnorm1, 
             params = list(DK_meanlog = param_defaults$disp_meanlog,
                           DK_sdlog = param_defaults$disp_sdlog,
                           SI_meanlog = param_defaults$serial_meanlog,
                           SI_sdlog = param_defaults$serial_sdlog), 
             cutoff = 0.95,
             N = 1000, 
             seed = 1345)
parallel::stopCluster(cl)

# Summarize the trees
links_all <- build_all_links(ttrees, N = 1000)
incs_all <- links_all[is.na(id_progen)]
links_consensus <- build_consensus_links(links_all)
tree_consensus <- build_consensus_tree(links_consensus, ttrees)

# Write out files of the trees and the links (consensus & all)
write_create(tree_consensus, here("analysis/out/ttrees/consensus_tree_ln0.95.csv"),
             fwrite)
write_create(links_consensus, here("analysis/out/ttrees/consensus_links_ln0.95.csv"),
             fwrite)
write_create(incs_all, here("analysis/out/ttrees/incs_all_ln0.95.csv"), 
       fwrite)

# write out incursions
load("data/sd_case_data.rda")
incs <- links_consensus[is.na(id_progen)]
incs %<>%
  as_tibble() %>%
  left_join(sd_case_data, by = c("id_case" = "id")) %>%
  mutate(date = dmy(symptoms_started)) %>%
  select(id = id_case, date,
         x_coord = utm_easting, y_coord = utm_northing, 
         prob)
write_create(incs, here("analysis/out/incursions.csv"), fwrite)

# Get priors
incs %>% 
  mutate(year = year(date)) %>% 
  group_by(year) %>% 
  summarize(n = n()) %$% mean(n)/52 -> incs_per_week
sd_case_data %>% 
  mutate(year = year(dmy(symptoms_started))) %>% 
  group_by(year) %>% summarize(n = n()) %$% mean(n)/52 -> cases_per_week

sd_incs <- cases_per_week/ 0.85 * 0.5 / 1.96 # 0.96

priors <- list(R0 = function(n) exp(rnorm(n, mean = 0.2, sd = 0.3)), # centered around 1.2
               iota = function(n) exp(rnorm(n, mean = 0, sd = 0.5)), # centered around 1
               k = function(n) exp(rnorm(n, mean = 0.25, sd = 0.25))) # centered around 1.35

write_create(priors, here("analysis/out/priors.rds"), saveRDS)
