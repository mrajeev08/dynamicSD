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
baseline <- fread("analysis/out/archive/campaign_sims/fixed_locs.csv")
limit_sspread <- fread("analysis/out/archive/int_sims/limit_sspreaders.csv")
limit_incs <- fread("analysis/out/archive/int_sims/limit_incs.csv")

summarize_ints <- function(ints) {
  
  ints %>%
    filter(!stopped) %>%
    select(max_conn_vill:mean_inc_total, cov, prop) %>%
    pivot_longer(max_conn_vill:mean_inc_total) %>%
    filter(cov != 0, prop != 0, 
           !(name %in% c("max_conn_vill", "mean_conn_vill", "nweeks_over_incs"))) %>%
    separate(name, into = c("type_stat", "class", "attr")) %>%
    mutate(stat = glue("{class} {attr}")) -> ints_summ
  
  out <- 
    ggplot(ints_summ, 
         aes(x = cov, y = value, color = prop, 
             group = interaction(cov, prop))) +
    geom_boxplot() +
    facet_wrap(stat ~ type_stat, scales = "free_y", ncol = 2) +
    ppars$theme_proj() +
    scale_color_distiller(palette = "BuPu") +
    labs(x = "Campaign coverage", 
         color = "Spatial coverage (proportion of \n villages vaccinated)")
  
  return(out)
}

plot_baseline_summ <- summarize_ints(baseline)
plot_sspread_summ <- summarize_ints(limit_sspread)
plot_incs_summ <- summarize_ints(limit_incs)

plot_conn_met <- 
  ggplot(filter(baseline, cov != 0, prop != 0, stopped == FALSE), 
         aes(x = cov, y = nweeks_over_incs, 
             color = max_conn_vill, 
             group = interaction(cov, prop))) +
  ggbeeswarm::geom_quasirandom() +
  facet_wrap(~prop) +
  labs(y = "Consecutive weeks above \n 95th% of introductions", 
       x = "Campaign coverage \n (proportion vaccinated in each village)", 
       color = "Peak connectivity of \n susceptible populations") +
  cowplot::theme_half_open(font_size = 12) +
  scale_color_distiller(palette = "BuPu") +
  labs()

# Interventions ----

# Do less than 2 months = probability of elimination!
baseline %>%
  group_by(cov, prop) %>%
  mutate(outbrk = nweeks_over_incs > 6 * 4 | stopped) %>%
  summarize(prop_outbrks = sum(outbrk)/n(),
            type = "Baseline") -> out

# limit transmission heterogeneity
limit_sspread %>%
  group_by(cov, prop) %>%
  mutate(outbrk = nweeks_over_incs > 6 * 4 | stopped) %>%
  summarize(prop_outbrks = sum(outbrk)/n(), 
            type = "Limit transmission heterogeneity") -> out_ss

# limit incursions
limit_incs %>%
  group_by(cov, prop) %>%
  mutate(outbrk = nweeks_over_incs > 6 * 4 | stopped) %>%
  summarize(prop_outbrks = sum(outbrk)/n(), 
            type = "Limit introductions") -> out_incs

comp_ints <- bind_rows(out, out_ss, out_incs)

plot_int_comparison <-
  ggplot(filter(comp_ints, cov != 0, prop != 0)) +
  geom_tile(aes(x = cov, y = prop, fill = prop_outbrks)) +
  facet_wrap(~type, labeller = labeller(type = label_wrap_gen(10))) +
  labs(x = "Campaign coverage", 
       y = "Spatial coverage (proportion of \n villages vaccinated)", 
       fill = "Outbreak \nprobability") +
  cowplot::theme_half_open(font_size = 12)

# hypothetical relationships ----
data.frame(prop_sus = seq(0, 1, by = 0.01)) %>%
  mutate(eff_prop = ifelse(prop_sus < 0.3, 0.3, prop_sus), 
         eff_red = 1 * eff_prop^2,
         type = "Limit introductions") -> inc_constraint

data.frame(prop_sus = seq(0, 1, by = 0.01)) %>%
  mutate(eff_prop = ifelse(prop_sus < 0.3, 0.3, prop_sus), 
         eff_red = 1 * eff_prop,
         type = "Limit transmission heterogeneity") -> het_constraint

constraints <- bind_rows(het_constraint, inc_constraint)

plot_constraint <-
  ggplot(constraints) +
  geom_line(aes(x = 1 - prop_sus, y = eff_red, color = type), 
            size = 2) + 
  ylim(c(0, 1)) +
  labs(x = "District vaccination coverage", y = "Effective reduction", 
       color = "Type of intervention") +
  scale_color_brewer(palette = "Accent") +
  ppars$theme_proj()


ggsave("analysis/figs/sfig_conn_mets.jpeg", plot_conn_met, height = 8, width = 8)
ggsave("analysis/figs/sfig_sspread_mets.jpeg", plot_sspread_summ)
ggsave("analysis/figs/sfig_incs_mets.jpeg", plot_incs_summ)
ggsave("analysis/figs/sfig_base_mets.jpeg", plot_baseline_summ)
ggsave("analysis/figs/mfig_int_comp.jpeg", plot_int_comparison)
ggsave("analysis/figs/sfig_conn_constraint.jpeg", plot_constraint)
