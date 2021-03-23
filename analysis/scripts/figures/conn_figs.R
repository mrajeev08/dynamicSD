# Connectivity sims 

library(data.table)
library(ggplot2)
library(foreach)
library(tibble)
library(readr)
library(dplyr)
library(magrittr)
library(glue)
library(ggtext)
library(tidytext)
library(tidyr)
source("R/utils.R")
source("R/figure_funs.R")

# params for plotting
ppars <- plotpars()

# Read in files -----
baseline <- fread("analysis/out/campaign_sims/fixed_locs.csv")
baseline <- fread("analysis/out/")

baseline %>%
  select(max_conn_vill:mean_inc_total, cov, prop) %>%
  pivot_longer(max_conn_vill:mean_inc_total) -> base_summ


ggplot(filter(base_summ, cov != 0, prop != 0), 
       aes(x = cov, y = value, color = prop, group = interaction(cov, prop))) +
 geom_boxplot() +
 facet_wrap(~name, scales = "free")

ggplot(filter(baseline, cov != 0, prop != 0), 
       aes(x = max_conn_vill, y = peak_chain_size, color = prop)) +
  geom_point() +
  facet_wrap(~cov)
ggplot(filter(baseline, cov != 0, prop != 0), 
       aes(x = max_conn_vill, y = nweeks_over_incs, color = prop)) +
  geom_point() +
  facet_wrap(~cov)

# Do less than 2 months = probability of elimination!


# Interventions ----

# hyp relationship

# simulated


